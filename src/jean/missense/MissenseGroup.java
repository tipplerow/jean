
package jean.missense;

import java.util.AbstractCollection;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import jam.lang.JamException;
import jam.util.CollectionUtil;
import jam.util.ListUtil;
import jam.util.ReadOnlyIterator;

import jean.ensembl.EnsemblTranscriptID;
import jean.hugo.HugoSymbol;
import jean.peptide.Peptide;
import jean.peptide.ProteinChange;
import jean.tcga.CellFraction;
import jean.tcga.TumorBarcode;

/**
 * Collects all missense mutations for a single tumor sample and gene.
 */
public final class MissenseGroup extends AbstractCollection<MissenseRecord> implements Comparable<MissenseGroup> {
    private final HugoSymbol hugoSymbol;
    private final TumorBarcode tumorBarcode;
    private final EnsemblTranscriptID transcriptID;

    private final Map<Integer, MissenseRecord> positionMap =
        new TreeMap<Integer, MissenseRecord>();

    private MissenseGroup(Collection<MissenseRecord> records) {
        MissenseRecord first = CollectionUtil.peek(records);

        if (first == null)
            throw JamException.runtime("No mutation records.");

        this.hugoSymbol = first.getHugoSymbol();
        this.tumorBarcode = first.getTumorBarcode();
        this.transcriptID = first.getTranscriptID();

        addRecords(records);
    }

    private void addRecords(Collection<MissenseRecord> records) {
        for (MissenseRecord record : records)
            addRecord(record);
    }

    private void addRecord(MissenseRecord record) {
        validateSymbol(record);
        validateBarcode(record);
        validateTranscript(record);

        MissenseRecord prev = putRecord(record);

        if (prev != null)
            putRecord(MissenseRecord.resolveDuplicate(prev, record));
    }

    private MissenseRecord putRecord(MissenseRecord record) {
        return positionMap.put(positionOf(record), record);        
    }

    private void validateSymbol(MissenseRecord record) {
        HugoSymbol groupSymbol = this.hugoSymbol;
        HugoSymbol recordSymbol = record.getHugoSymbol();

        if (!recordSymbol.equals(groupSymbol))
            throw JamException.runtime("Inconsistent HUGO symbols: [%s != %s].", recordSymbol, groupSymbol);
    }

    private void validateBarcode(MissenseRecord record) {
        TumorBarcode groupBarcode = this.tumorBarcode;
        TumorBarcode recordBarcode = record.getTumorBarcode();

        if (!recordBarcode.equals(groupBarcode))
            throw JamException.runtime("Inconsistent barcodes: [%s != %s].", recordBarcode, groupBarcode);
    }

    private void validateTranscript(MissenseRecord record) {
        EnsemblTranscriptID groupTranscript = this.transcriptID;
        EnsemblTranscriptID recordTranscript = record.getTranscriptID();
        
        if (groupTranscript == null && recordTranscript != null)
            throw JamException.runtime("Mixed null and non-null transcripts.");

        if (groupTranscript != null && recordTranscript == null)
            throw JamException.runtime("Mixed null and non-null transcripts.");

        if (groupTranscript != null && !recordTranscript.equals(groupTranscript))
            throw JamException.runtime("Inconsistent transcripts: [%s != %s].", recordTranscript, groupTranscript);
    }

    private static int positionOf(MissenseRecord record) {
        return record.getProteinChange().getPosition();
    }

    /**
     * A comparator that orders missense groups by tumor barcode
     * first, then HUGO symbol.
     */
    public static Comparator<MissenseGroup> COMPARATOR =
        new Comparator<MissenseGroup>() {
            @Override public int compare(MissenseGroup group1, MissenseGroup group2) {
                int cmp = group1.tumorBarcode.compareTo(group2.tumorBarcode);
                
                if (cmp != 0)
                    return cmp;
                else
                    return group1.hugoSymbol.compareTo(group2.hugoSymbol);
            }
        };

    /**
     * Creates a new group with a fixed collection of mutations.
     *
     * @param records the missense mutations for a single tumor sample
     * and gene.
     *
     * @return a new group with the specified missense mutations.
     *
     * @throws RuntimeException unless all mutation records refer to
     * the same tumor sample and gene and any conflicting duplicate
     * mutations can be resolved (by keeping the one with the greatest
     * cancer cell fraction).
     */
    public static MissenseGroup create(Collection<MissenseRecord> records) {
        return new MissenseGroup(records);
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
     * Returns a list of the protein changes in this group.
     *
     * @return a list of the protein changes in this group.
     */
    public List<ProteinChange> getProteinChanges() {
        return ListUtil.apply(positionMap.values(), record -> record.getProteinChange());
    }

    /**
     * Applies the mutations in this group to the native peptide.
     *
     * @param native_ the native peptide.
     *
     * @return the mutated peptide.
     *
     * @throws RuntimeException unless the native peptide is valid
     * (the residues at each mutation position in the native peptide
     * match the native residues in the peptide change data).
     */
    public Peptide mutate(Peptide native_) {
        return native_.mutate(getProteinChanges());
    }

    /**
     * Returns a read-only view of the (unit-offset) locations of the
     * mutations in this group.
     *
     * @return a read-only view of the (unit-offset) locations of the
     * mutations in this group.
     */
    public Set<Integer> viewPositions() {
        return Collections.unmodifiableSet(positionMap.keySet());
    }

    /**
     * Compares this group to another by tumor barcode first, then
     * HUGO symbol.
     *
     * @return an integer less than, equal to, or greater than one
     * according to whether this group should be ordered before, in
     * place, or after the input group.
     */
    @Override public int compareTo(MissenseGroup that) {
        return COMPARATOR.compare(this, that);
    }

    /**
     * Returns an iterator over the records in this group.
     *
     * @return an iterator over the records in this group.
     */
    @Override public Iterator<MissenseRecord> iterator() {
            return ReadOnlyIterator.create(positionMap.values());
    }

    /**
     * Returns the number of records in this group.
     *
     * @return the number of records in this group.
     */
    @Override public int size() {
        return positionMap.size();
    }

    @Override public String toString() {
        return positionMap.values().toString();
    }
}
