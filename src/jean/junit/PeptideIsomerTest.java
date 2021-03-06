
package jean.junit;

import jean.peptide.Peptide;
import jean.peptide.PeptideIsomer;
import jean.peptide.Residue;

import org.junit.*;
import static org.junit.Assert.*;

public class PeptideIsomerTest {
    @Test public void testIsomerKey() {
        Peptide pep = Peptide.of(Residue.His,
                                 Residue.Cys,
                                 Residue.Gln,
                                 Residue.Gln,
                                 Residue.Cys,
                                 Residue.Gln,
                                 Residue.Ala,
                                 Residue.Lys);

        assertEquals("ACCHKQQQ", PeptideIsomer.isomerKey(pep));
    }

    @Test public void testMapIsomers() {
        assertEquals( 210, PeptideIsomer.mapIsomers(2).elementSet().size());
        assertEquals(1540, PeptideIsomer.mapIsomers(3).elementSet().size());
        assertEquals(8855, PeptideIsomer.mapIsomers(4).elementSet().size());
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("jean.junit.PeptideIsomerTest");
    }
}
