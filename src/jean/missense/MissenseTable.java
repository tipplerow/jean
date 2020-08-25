
package jean.missense;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import jam.app.JamLogger;
import jam.util.PairKeyTable;

import jean.hugo.HugoSymbol;
import jean.tcga.TumorBarcode;

/**
 * Indexes missesnse mutations by tumor barcode and HUGO symbol.
 */
public final class MissenseTable {
    //
    // Class-specific containers...
    //
    private static final class RecordList extends ArrayList<MissenseRecord> {}

    // All mutations indexed by barcode (outer) and symbol (inner)...
    private final PairKeyTable<TumorBarcode, HugoSymbol, RecordList> table = PairKeyTable.hash();

    private MissenseTable(Collection<MissenseRecord> records) {
        fillMap(records);
    }

    private void fillMap(Collection<MissenseRecord> records) {
        for (MissenseRecord record : records)
            addRecord(record);
    }

    private void addRecord(MissenseRecord record) {
        HugoSymbol   symbol  = record.getHugoSymbol();
        TumorBarcode barcode = record.getTumorBarcode();

        recordList(barcode, symbol).add(record);
    }

    private RecordList recordList(TumorBarcode barcode, HugoSymbol symbol) {
        RecordList recordList = table.get(barcode, symbol);

        if (recordList == null) {
            recordList = new RecordList();
            table.put(barcode, symbol, recordList);
        }

        return recordList;
    }

    /**
     * Populates a table by reading all missense mutation records from
     * a given file.
     *
     * @param fileName the path to the missense mutation file.
     *
     * @return a table containing all missense mutation records in the
     * given file.
     *
     * @throws RuntimeException unless the file can be opened for
     * reading and contains properly formatted records.
     */
    public static MissenseTable load(String fileName) {
        List<MissenseRecord> records = MissenseParser.parse(fileName);
        JamLogger.info("MissenseTable: Loaded [%d] records.", records.size());
        
        return load(records);
    }

    /**
     * Populates a table from a collection of missense mutation
     * records.
     *
     * @param records the records to be indexed in the table.
     *
     * @return a table containing all missense mutation records in
     * the given collection.
     *
     * @throws RuntimeException unless all Ensembl transcripts are
     * identical for records with the same tumor sample and gene.
     */
    public static MissenseTable load(Collection<MissenseRecord> records) {
        return new MissenseTable(records);
    }

    /**
     * Identifies tumor-gene pairs contained in this mutation table.
     *
     * @param barcode the tumor barcode of interest.
     *
     * @param symbol the gene of interest.
     *
     * @return {@code true} iff this table contains mutations for the
     * specified tumor-gene pair.
     */
    public boolean contains(TumorBarcode barcode, HugoSymbol symbol) {
        return table.contains(barcode, symbol);
    }

    /**
     * Counts the total number of missense mutations in a given tumor.
     *
     * @param barcode the tumor barcode of interest.
     *
     * @return the total number of missense mutations in the specified
     * tumor.
     */
    public int count(TumorBarcode barcode) {
        int total = 0;

        for (HugoSymbol symbol : table.viewInnerKeys(barcode))
            total += count(barcode, symbol);

        return total;
    }

    /**
     * Counts the number of missense mutations for a given tumor
     * and gene.
     *
     * @param barcode the tumor barcode of interest.
     *
     * @param symbol the gene of interest.
     *
     * @return the number of missense mutations for the specified
     * tumor and gene.
     */
    public int count(TumorBarcode barcode, HugoSymbol symbol) {
        return lookup(barcode, symbol).size();
    }

    /**
     * Returns all missense mutations for a given tumor and gene.
     *
     * @param barcode the tumor barcode of interest.
     *
     * @param symbol the HUGO symbol of interest.
     *
     * @return an immutable list containing all missense mutations for
     * the specified tumor and gene (or an empty list if there are no
     * matching mutations).
     */
    public List<MissenseRecord> lookup(TumorBarcode barcode, HugoSymbol symbol) {
        RecordList recordList = table.get(barcode, symbol);

        if (recordList != null)
            return Collections.unmodifiableList(recordList);
        else
            return Collections.emptyList();
    }

    /**
     * Returns a read-only view of all tumor barcodes in this table.
     *
     * @return a read-only view of all tumor barcodes in this table.
     */
    public Set<TumorBarcode> viewBarcodes() {
        return table.viewOuterKeys();
    }

    /**
     * Returns a read-only view of all mutated genes for a given tumor.
     *
     * @param barcode the tumor barcode of interest.
     *
     * @return a read-only view of all mutated genes for the specified
     * tumor.
     */
    public Set<HugoSymbol> viewSymbols(TumorBarcode barcode) {
        return table.viewInnerKeys(barcode);
    }
}
