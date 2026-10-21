/*
 * Multi-threaded drop-in replacement for Picard CollectRnaSeqMetrics.
 *
 * Accepts exactly the same command line arguments as `picard CollectRnaSeqMetrics`.
 * The number of threads is given by the system property `threads` (default: 1).
 *
 *   java -Dthreads=4 -cp picard.jar ParallelRnaSeqMetrics.java --REF_FLAT ... -I in.bam -O out.txt
 *
 * With threads > 1, the coordinate-sorted and indexed input BAM is split into
 * genomic chunks (balanced by compressed BAM size, using the BAM index). Every
 * record is assigned to exactly one chunk by its alignment start, and the
 * unplaced unmapped reads form an extra chunk. Each chunk is processed by
 * Picard's own RnaSeqMetricsCollector, in parallel. At the end the per-chunk
 * results are merged: counts are summed and per-transcript coverage arrays
 * are added element-wise. The coverage based metrics (MEDIAN_CV_COVERAGE,
 * 5'/3' biases, normalized coverage histogram) are then computed by Picard on
 * the merged data, so the output is the same as the single-threaded run.
 *
 * Falls back to the standard single-threaded Picard code when threads <= 1,
 * or when the input is not an indexed BAM, or with STOP_AFTER, REFERENCE_SEQUENCE
 * or a METRIC_ACCUMULATION_LEVEL other than ALL_READS.
 */

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.samtools.util.BlockGunzipper;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.analysis.CollectRnaSeqMetrics;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.RnaSeqMetrics;
import picard.analysis.directed.RnaSeqMetricsCollector;
import picard.annotation.Gene;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.metrics.PerUnitMetricCollector;
import com.intel.gkl.compression.IntelInflaterFactory;

import java.io.File;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

@CommandLineProgramProperties(
        summary = "Multi-threaded CollectRnaSeqMetrics (set the number of threads with -Dthreads=N)",
        oneLineSummary = "Multi-threaded CollectRnaSeqMetrics",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
public class ParallelRnaSeqMetrics extends CollectRnaSeqMetrics {

    private static final Log LOG = Log.getInstance(ParallelRnaSeqMetrics.class);

    /** Number of chunks per thread, small chunks give a better load balancing. */
    private static final int CHUNKS_PER_THREAD = 8;
    /** Do not split regions below this size (the BAM linear index resolution is 16 kb). */
    private static final int MIN_REGION_LENGTH = 1 << 15;

    /** The counts of RnaSeqMetrics, all other fields are derived from these by Picard. */
    private static final String[] COUNT_FIELDS = {
            "PF_BASES", "PF_ALIGNED_BASES", "CODING_BASES", "UTR_BASES", "INTRONIC_BASES", "INTERGENIC_BASES",
            "IGNORED_READS", "CORRECT_STRAND_READS", "INCORRECT_STRAND_READS",
            "NUM_R1_TRANSCRIPT_STRAND_READS", "NUM_R2_TRANSCRIPT_STRAND_READS", "NUM_UNEXPLAINED_READS"};

    /** A chunk of the input: all records with alignment start in [start, end] on contig refIndex,
     *  or the unplaced unmapped reads if refIndex == -1. */
    private record Region(int refIndex, String contig, int start, int end) {}

    private record ChunkResult(Region region, RnaSeqMetricsCollector collector, long records) {}

    public static void main(final String[] argv) {
        final int threads = Integer.getInteger("threads", 1);
        System.exit(threads <= 1
                ? new CollectRnaSeqMetrics().instanceMain(argv)
                : new ParallelRnaSeqMetrics().instanceMainParallel(argv, threads));
    }

    /** Same as CommandLineProgram.instanceMain(), but runs the parallel code instead of doWork(). */
    private int instanceMainParallel(final String[] argv, final int threads) {
        if (!parseArgs(argv)) return 1;

        // Same default headers as Picard CollectRnaSeqMetrics writes into the metrics file
        final Date startDate = new Date();
        final String commandLine = getCommandLine().replaceFirst("^" + getClass().getSimpleName(),
                CollectRnaSeqMetrics.class.getSimpleName());
        getDefaultHeaders().add(new StringHeader(commandLine));
        getDefaultHeaders().add(new StringHeader("Started on: " + startDate));

        Log.setGlobalLogLevel(VERBOSITY);
        if (!USE_JDK_INFLATER) BlockGunzipper.setDefaultInflaterFactory(new IntelInflaterFactory());
        SamReaderFactory.setDefaultValidationStringency(VALIDATION_STRINGENCY);
        if (!QUIET) System.err.println("[" + startDate + "] " + commandLine + " (threads: " + threads + ")");

        if (OUTPUT == null) OUTPUT = output.getOutputFile();

        final int ret;
        if (STOP_AFTER > 0 || REFERENCE_SEQUENCE != null
                || !METRIC_ACCUMULATION_LEVEL.equals(Collections.singleton(MetricAccumulationLevel.ALL_READS))) {
            LOG.warn("STOP_AFTER, REFERENCE_SEQUENCE or METRIC_ACCUMULATION_LEVEL other than ALL_READS: running single-threaded");
            ret = doWork();
        } else {
            ret = runParallel(threads);
        }

        if (!QUIET) {
            final Date endDate = new Date();
            System.err.printf("[%s] %s done. Elapsed time: %.2f minutes.%n", endDate, getClass().getName(),
                    (endDate.getTime() - startDate.getTime()) / 60000d);
        }
        return ret;
    }

    private int runParallel(final int threads) {
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);
        final SAMFileHeader header = in.getFileHeader();
        if (!in.hasIndex() || in.indexing().getIndex() == null) {
            CloserUtil.close(in);
            LOG.warn("Input is not an indexed BAM file: running single-threaded");
            return doWork();
        }
        if (header.getSortOrder() != SAMFileHeader.SortOrder.coordinate && !ASSUME_SORTED) {
            throw new PicardException("File " + INPUT.getAbsolutePath() + " should be coordinate sorted but " +
                    "the header says the sort order is " + header.getSortOrder() + ". If you believe the file " +
                    "to be coordinate sorted you may pass ASSUME_SORTED=true");
        }

        // Builds Picard's collector (gene/rRNA overlap detectors) exactly as CollectRnaSeqMetrics does
        setup(header, INPUT);
        final RnaSeqMetricsCollector mainCollector = getField(this, CollectRnaSeqMetrics.class, "collector");

        final BAMIndex index = in.indexing().getIndex();
        final List<Region> regions = makeRegions(header, index, threads);
        final long[] indexedCounts = new long[header.getSequenceDictionary().size()];
        for (int i = 0; i < indexedCounts.length; i++) {
            final BAMIndexMetaData md = index.getMetaData(i);
            indexedCounts[i] = md == null ? 0 : (long) md.getAlignedRecordCount() + md.getUnalignedRecordCount();
        }
        CloserUtil.close(in);
        LOG.info("Processing " + regions.size() + " chunks with " + threads + " threads");

        // One reader per thread, one collector per chunk (gene and rRNA overlap detectors are shared, read-only)
        final List<SamReader> readers = Collections.synchronizedList(new ArrayList<>());
        final ThreadLocal<SamReader> reader = ThreadLocal.withInitial(() -> {
            final SamReader r = SamReaderFactory.makeDefault().open(INPUT);
            readers.add(r);
            return r;
        });
        final ExecutorService pool = Executors.newFixedThreadPool(threads);
        final List<Future<ChunkResult>> futures = new ArrayList<>();
        for (final Region region : regions) {
            futures.add(pool.submit(() -> processRegion(reader.get(), region, newCollector(mainCollector, header))));
        }
        pool.shutdown();

        // Merge the chunks in genomic order, while the remaining ones are still running
        final PerUnitMetricCollector<RnaSeqMetrics, Integer, SAMRecord> merged = mainCollector.getAllReadsCollector();
        final long[] processedCounts = new long[indexedCounts.length];
        long total = 0;
        try {
            for (int i = 0; i < futures.size(); i++) {
                final ChunkResult chunk = futures.get(i).get();
                futures.set(i, null);
                if (chunk.region.refIndex >= 0) processedCounts[chunk.region.refIndex] += chunk.records;
                total += chunk.records;
                mergeInto(merged, chunk.collector.getAllReadsCollector());
            }
        } catch (final Exception e) {
            pool.shutdownNow();
            throw new PicardException("Failed to process chunk: " + e.getMessage(), e);
        } finally {
            readers.forEach(CloserUtil::close);
        }

        // Sanity check: every placed record has been processed exactly once
        if (Arrays.stream(indexedCounts).sum() > 0) {
            for (int i = 0; i < indexedCounts.length; i++) {
                if (indexedCounts[i] != processedCounts[i]) {
                    throw new PicardException("Record count mismatch for " + header.getSequence(i).getSequenceName() +
                            ": index " + indexedCounts[i] + ", processed " + processedCounts[i]);
                }
            }
        } else {
            LOG.warn("BAM index has no record counts, skipping record count check");
        }
        LOG.info("Processed " + total + " records.");

        // Picard computes the derived and coverage based metrics and writes the metrics file
        finish();
        return 0;
    }

    /** Splits each contig into regions with roughly equal compressed BAM size, plus one region for unplaced reads. */
    private static List<Region> makeRegions(final SAMFileHeader header, final BAMIndex index, final int threads) {
        final List<SAMSequenceRecord> sequences = header.getSequenceDictionary().getSequences();
        long totalBytes = 0;
        for (final SAMSequenceRecord seq : sequences) {
            totalBytes += spanBytes(index, seq.getSequenceIndex(), 1, seq.getSequenceLength());
        }
        final long targetBytes = Math.max(1, totalBytes / ((long) threads * CHUNKS_PER_THREAD));

        final List<Region> regions = new ArrayList<>();
        for (final SAMSequenceRecord seq : sequences) {
            split(index, seq, 1, seq.getSequenceLength(), targetBytes, regions);
        }
        regions.add(new Region(-1, null, 0, 0));
        return regions;
    }

    private static void split(final BAMIndex index, final SAMSequenceRecord seq, final int start, final int end,
                              final long targetBytes, final List<Region> regions) {
        final long bytes = spanBytes(index, seq.getSequenceIndex(), start, end);
        if (bytes == 0) return;
        if (bytes <= targetBytes || end - start + 1 <= MIN_REGION_LENGTH) {
            regions.add(new Region(seq.getSequenceIndex(), seq.getSequenceName(), start, end));
            return;
        }
        final int mid = start + (end - start) / 2;
        split(index, seq, start, mid, targetBytes, regions);
        split(index, seq, mid + 1, end, targetBytes, regions);
    }

    /** Approximate compressed size of the records overlapping a region, from the BAM index. */
    private static long spanBytes(final BAMIndex index, final int refIndex, final int start, final int end) {
        long bytes = 0;
        for (final Chunk chunk : index.getSpanOverlapping(refIndex, start, end).getChunks()) {
            bytes += (chunk.getChunkEnd() >>> 16) - (chunk.getChunkStart() >>> 16) + 1;
        }
        return bytes;
    }

    private static ChunkResult processRegion(final SamReader reader, final Region region,
                                             final RnaSeqMetricsCollector collector) {
        long records = 0;
        try (SAMRecordIterator it = region.refIndex < 0
                ? reader.queryUnmapped()
                : reader.queryOverlapping(region.contig, region.start, region.end)) {
            while (it.hasNext()) {
                final SAMRecord rec = it.next();
                // Assign each record to exactly one region: the one containing its alignment start
                if (region.refIndex >= 0 && (rec.getAlignmentStart() < region.start || rec.getAlignmentStart() > region.end)) {
                    continue;
                }
                collector.acceptRecord(rec, null);
                records++;
            }
        }
        return new ChunkResult(region, collector, records);
    }

    /** A new collector with the same settings and the same (shared) overlap detectors as the main collector. */
    private RnaSeqMetricsCollector newCollector(final RnaSeqMetricsCollector main, final SAMFileHeader header) {
        final OverlapDetector<Gene> genes = getField(main, RnaSeqMetricsCollector.class, "geneOverlapDetector");
        final OverlapDetector<htsjdk.samtools.util.Interval> ribosomal =
                getField(main, RnaSeqMetricsCollector.class, "ribosomalSequenceOverlapDetector");
        final Set<Integer> ignored = getField(main, RnaSeqMetricsCollector.class, "ignoredSequenceIndices");
        return new RnaSeqMetricsCollector(METRIC_ACCUMULATION_LEVEL, header.getReadGroups(),
                RIBOSOMAL_INTERVALS != null ? 0L : null, genes, ribosomal, new HashSet<>(ignored),
                MINIMUM_LENGTH, STRAND_SPECIFICITY, RRNA_FRAGMENT_PERCENTAGE, true, END_BIAS_BASES);
    }

    /** Adds the counts and the per-transcript coverage of one chunk to the merged collector. */
    private static void mergeInto(final PerUnitMetricCollector<RnaSeqMetrics, Integer, SAMRecord> target,
                                  final PerUnitMetricCollector<RnaSeqMetrics, Integer, SAMRecord> source) {
        final Class<?> perUnit = target.getClass();
        final RnaSeqMetrics t = getField(target, perUnit, "metrics");
        final RnaSeqMetrics s = getField(source, perUnit, "metrics");
        try {
            for (final String name : COUNT_FIELDS) {
                final Field f = RnaSeqMetrics.class.getField(name);
                f.setLong(t, f.getLong(t) + f.getLong(s));
            }
        } catch (final ReflectiveOperationException e) {
            throw new PicardException("Cannot merge metrics", e);
        }
        if (t.RIBOSOMAL_BASES != null) t.RIBOSOMAL_BASES += s.RIBOSOMAL_BASES;

        final Map<Gene.Transcript, int[]> tCov = getField(target, perUnit, "coverageByTranscript");
        final Map<Gene.Transcript, int[]> sCov = getField(source, perUnit, "coverageByTranscript");
        for (final Map.Entry<Gene.Transcript, int[]> e : sCov.entrySet()) {
            final int[] cov = tCov.get(e.getKey());
            if (cov == null) {
                tCov.put(e.getKey(), e.getValue());
            } else {
                final int[] add = e.getValue();
                for (int i = 0; i < cov.length; i++) cov[i] += add[i];
            }
        }
    }

    @SuppressWarnings("unchecked")
    private static <T> T getField(final Object obj, final Class<?> cls, final String name) {
        try {
            final Field f = cls.getDeclaredField(name);
            f.setAccessible(true);
            return (T) f.get(obj);
        } catch (final ReflectiveOperationException e) {
            throw new PicardException("Unsupported Picard version, cannot access " + cls.getName() + "." + name, e);
        }
    }
}
