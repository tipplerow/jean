
package jean.junit;

import java.util.List;

import jam.math.IntRange;

import jean.ensembl.EnsemblProteinDb;
import jean.hugo.HugoMaster;
import jean.missense.MissenseGroup;
import jean.missense.MissenseTable;
import jean.neo.PeptidePairEngine;
import jean.neo.PeptidePairRecord;
import jean.peptide.Peptide;
import jean.tcga.TumorBarcode;

import org.junit.*;
import static org.junit.Assert.*;

public class PeptidePairEngineTest {
    @Test public void testGenerate() {
        HugoMaster hugoMaster = HugoMaster.load("data/test/hugo_master_test.tsv");
        EnsemblProteinDb ensemblDb = EnsemblProteinDb.load("data/test/ensembl_test2.fa");

        PeptidePairEngine.initialize(hugoMaster, ensemblDb);

        TumorBarcode barcode1 = TumorBarcode.instance("barcode1");
        TumorBarcode barcode2 = TumorBarcode.instance("barcode2");

        MissenseTable missenseTable = MissenseTable.load("data/test/ppe_missense.maf");
        List<MissenseGroup> missenseGroups = missenseTable.group(barcode1);

        assertEquals(1, missenseGroups.size());

        List<PeptidePairRecord> pairRecords =
            PeptidePairEngine.generate(missenseGroups.get(0), 9);

        assertEquals(21, pairRecords.size());

        assertEquals(IntRange.instance(9, 17), pairRecords.get(0).getPeptideRange());
        assertEquals(IntRange.instance(20, 28), pairRecords.get(11).getPeptideRange());
        assertEquals(IntRange.instance(168, 176), pairRecords.get(12).getPeptideRange());
        assertEquals(IntRange.instance(176, 184), pairRecords.get(20).getPeptideRange());

        missenseGroups = missenseTable.group(barcode2);
        pairRecords = PeptidePairEngine.generate(missenseGroups.get(0), 9);

        assertEquals(1, missenseGroups.size());
        assertEquals(16, pairRecords.size());

        assertEquals(IntRange.instance(1, 9), pairRecords.get(0).getPeptideRange());
        assertEquals(IntRange.instance(3, 11), pairRecords.get(2).getPeptideRange());
        assertEquals(IntRange.instance(168, 176), pairRecords.get(3).getPeptideRange());
        assertEquals(IntRange.instance(181, 189), pairRecords.get(15).getPeptideRange());

        pairRecords = PeptidePairEngine.generate(missenseTable, 9);

        for (PeptidePairRecord record : pairRecords)
            System.out.println(record);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("jean.junit.PeptidePairEngineTest");
    }
}
