
package jean.missense;

import java.util.List;
import java.util.Set;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

import jam.util.CollectionUtil;
import jam.util.ListUtil;
import jam.util.MultisetUtil;

import jean.ensembl.EnsemblTranscriptID;
import jean.hugo.HugoSymbol;
import jean.peptide.ProteinChange;
import jean.tcga.TumorBarcode;
import jean.tcga.CellFraction;

/**
 * Encapsulates information contained in one record from a Mutation
 * Annotation Format (MAF) file that describes only <b>missense</b>
 * mutations.
 */
public final class MissenseRecord {
    private final TumorBarcode        tumorBarcode;
    private final EnsemblTranscriptID transcriptID;
    private final HugoSymbol          hugoSymbol;
    private final ProteinChange       proteinChange;
    private final CellFraction        cellFraction;

    /**
     * Creates a new missense mutation record.
     *
     * @param tumorBarcode the tumor in which the mutation occurred.
     *
     * @param transcriptID the Ensembl identifier for the mutated RNA
     * transcript.
     *
     * @param hugoSymbol the HUGO symbol for the mutated gene.
     *
     * @param proteinChange the description of the single-residue change.
     *
     * @param cellFraction the cancer cell fraction of the mutation.
     */
    public MissenseRecord(TumorBarcode        tumorBarcode,
                          EnsemblTranscriptID transcriptID,
                          HugoSymbol          hugoSymbol,
                          ProteinChange       proteinChange,
                          CellFraction        cellFraction) {
        this.tumorBarcode  = tumorBarcode;
        this.transcriptID  = transcriptID;
        this.hugoSymbol    = hugoSymbol;
        this.proteinChange = proteinChange;
        this.cellFraction  = cellFraction;
    }

    /**
     * Filters a list of missense records and retains only those with
     * cell fractions above a threshold.
     *
     * @param records the records to filter.
     *
     * @param ccfThreshold the cancer cell fraction threshold.
     *
     * @return a list containing only records with cell fractions
     * above the specified threshold.
     */
    public static List<MissenseRecord> filterCellFraction(List<MissenseRecord> records, CellFraction ccfThreshold) {
        return ListUtil.filter(records, record -> record.getCellFraction().above(ccfThreshold));
    }

    /**
     * Filters a list of missense records and retains only those with
     * Ensembl transcripts that match the most common transcript.
     *
     * <p>If all records have the same transcript, then the input list
     * is returned unaltered.  If the list does not contain a unique
     * primary transcript (no single transcript occurs more than any
     * other), then an empty list is returned.
     *
     * @param records the records to filter.
     *
     * @return a list containing only those records having the primary
     * transcript.
     */
    public static List<MissenseRecord> filterPrimaryTranscript(List<MissenseRecord> records) {
        if (records.size() < 2)
            return records;

        Multiset<EnsemblTranscriptID> transcripts = HashMultiset.create();

        for (MissenseRecord record : records)
            transcripts.add(record.getTranscriptID());

        if (MultisetUtil.countUnique(transcripts) == 1)
            return records;

        Set<EnsemblTranscriptID> primary = MultisetUtil.mostCommon(transcripts);

        if (primary.size() == 1)
            return filterTranscript(records, CollectionUtil.peek(primary));
        else
            return List.of();
    }

    /**
     * Filters a list of missense records and retains only those whose
     * Ensembl transcript matches a target.
     *
     * @param records the records to filter.
     *
     * @param transcriptID the target transcript to match.
     *
     * @return a list containing only those records whose transcript
     * matches the target.
     */
    public static List<MissenseRecord> filterTranscript(List<MissenseRecord> records, EnsemblTranscriptID transcriptID) {
        if (transcriptID == null)
            return ListUtil.filter(records, record -> !record.hasTranscriptID());
        else
            return ListUtil.filter(records, record -> record.getTranscriptID().equals(transcriptID));
    }

    /**
     * Returns the HUGO symbol of the mutated gene.
     *
     * @return the HUGO symbol of the mutated gene.
     */
    public HugoSymbol getHugoSymbol() {
        return hugoSymbol;
    }

    /**
     * Returns the description of the single-residue change.
     *
     * @return the description of the single-residue change.
     */
    public ProteinChange getProteinChange() {
        return proteinChange;
    }

    /**
     * Returns the Ensembl identifier for the mutated RNA transcript.
     *
     * @return the Ensembl identifier for the mutated RNA transcript.
     */
    public EnsemblTranscriptID getTranscriptID() {
        return transcriptID;
    }

    /**
     * Returns the tumor in which the mutation occurred.
     *
     * @return the tumor in which the mutation occurred.
     */
    public TumorBarcode getTumorBarcode() {
        return tumorBarcode;
    }

    /**
     * Returns the cancer cell fraction for the mutation.
     *
     * @return the cancer cell fraction for the mutation.
     */
    public CellFraction getCellFraction() {
        return cellFraction;
    }

    /**
     * Identifies records with non-{@code null} Ensembl transcript
     * identifiers.
     *
     * @return {@code true} iff this record has a non-{@code null}
     * Ensembl transcript identifier.
     */
    public boolean hasTranscriptID() {
        return transcriptID != null;
    }

}
