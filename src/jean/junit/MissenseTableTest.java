
package jean.junit;

import jean.hugo.HugoSymbol;
import jean.missense.MissenseRecord;
import jean.missense.MissenseTable;
import jean.tcga.TumorBarcode;

import org.junit.*;
import static org.junit.Assert.*;

public class MissenseTableTest {
    private static final String TCGA_MAF = "data/test/TCGA_Missense.maf";

    private static final TumorBarcode barcode1 = TumorBarcode.instance("TCGA-02-0003-01A");
    private static final TumorBarcode barcode2 = TumorBarcode.instance("TCGA-02-0033-01A");
    private static final TumorBarcode barcode3 = TumorBarcode.instance("TCGA-02-0047-01A");

    private static final HugoSymbol ABR     = HugoSymbol.instance("ABR");
    private static final HugoSymbol ACADS   = HugoSymbol.instance("ACADS");
    private static final HugoSymbol ADAMTS2 = HugoSymbol.instance("ADAMTS2");
    private static final HugoSymbol GPR158  = HugoSymbol.instance("GPR158");
    private static final HugoSymbol OR5M3 = HugoSymbol.instance("OR5M3");
    private static final HugoSymbol ZNF583  = HugoSymbol.instance("ZNF583");

    @Test public void testTCGA() {
        MissenseTable table = MissenseTable.load(TCGA_MAF);

        assertTrue(table.viewBarcodes().contains(barcode1));
        assertTrue(table.viewBarcodes().contains(barcode2));
        assertTrue(table.viewBarcodes().contains(barcode3));
        assertFalse(table.viewBarcodes().contains(TumorBarcode.instance("no such")));

        assertEquals(44, table.count(barcode1));
        assertEquals(23, table.count(barcode2));
        assertEquals(33, table.count(barcode3));
        assertEquals(0,  table.count(TumorBarcode.instance("no such")));

        assertTrue(table.contains(barcode1, ZNF583));
        assertTrue(table.contains(barcode2, ACADS));
        assertFalse(table.contains(barcode1, ACADS));
        assertFalse(table.contains(barcode2, ZNF583));

        // Not 33 because there are two genes with double mutations...
        assertEquals(31, table.viewSymbols(barcode3).size());

        assertEquals(1, table.lookup(barcode3, ABR).size());
        assertEquals(2, table.lookup(barcode3, ADAMTS2).size());
        assertEquals(2, table.lookup(barcode3, GPR158).size());
        assertEquals(1, table.lookup(barcode3, OR5M3).size());

        assertEquals(0, table.lookup(barcode2, OR5M3).size());
        assertEquals(0, table.lookup(barcode3, ACADS).size());

        assertRecord(table.lookup(barcode3, ABR).get(0),     barcode3, ABR,     "ENST00000302538", "G532S");
        assertRecord(table.lookup(barcode3, ADAMTS2).get(0), barcode3, ADAMTS2, "ENST00000251582", "T805M");
        assertRecord(table.lookup(barcode3, ADAMTS2).get(1), barcode3, ADAMTS2, "ENST00000251582", "D361N");
        assertRecord(table.lookup(barcode3, GPR158).get(0),  barcode3, GPR158,  "ENST00000376351", "D778Y");
        assertRecord(table.lookup(barcode3, GPR158).get(1),  barcode3, GPR158,  "ENST00000376351", "G784R");
        assertRecord(table.lookup(barcode3, OR5M3).get(0),   barcode3, OR5M3,   "ENST00000312240", "T153M");
    }

    private void assertRecord(MissenseRecord record,
                              TumorBarcode   tumorBarcode,
                              HugoSymbol     hugoSymbol,
                              String         transcriptID,
                              String         proteinChange) {
        assertEquals(tumorBarcode,  record.getTumorBarcode());
        assertEquals(hugoSymbol,    record.getHugoSymbol());
        assertEquals(transcriptID,  record.getTranscriptID().getKey());
        assertEquals(proteinChange, record.getProteinChange().format());
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("jean.junit.MissenseTableTest");
    }
}
