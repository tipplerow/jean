
package jean.junit;

import jean.neo.NeoPeptide;
import jean.neo.PeptidePair;
import jean.neo.PeptideType;
import jean.neo.SelfPeptide;
import jean.peptide.Peptide;

import org.junit.*;
import static org.junit.Assert.*;

public class PeptidePairTest {
    @Test public void testInstance() {
        Peptide p1 = Peptide.instance("QVSRDQVLD");
        Peptide p2 = Peptide.instance("QVSRDQVLE");

        SelfPeptide self = SelfPeptide.instance(p1);
        NeoPeptide  neo  = NeoPeptide.instance(p2);
        PeptidePair pair = PeptidePair.instance(self, neo);

        assertEquals(p1, pair.self());
        assertEquals(p2, pair.neo());

        assertEquals(PeptideType.NEO, neo.getType());
        assertEquals(PeptideType.SELF, self.getType());
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("jean.junit.PeptidePairTest");
    }
}
