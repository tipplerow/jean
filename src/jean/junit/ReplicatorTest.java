
package jean.junit;

import java.util.ArrayList;
import java.util.List;

import jean.bio.Replicator;

import org.junit.*;
import static org.junit.Assert.*;

public class ReplicatorTest {
    private static final Replicator FIRST = Replicator.create();

    @Test public void testFirst() {
        assertEquals(0, FIRST.getIndex());
        assertEquals(null, FIRST.getParent());
        assertTrue(FIRST.isFounder());
    }

    @Test public void testDaughter() {
        Replicator step1 = FIRST.replicate();

        assertEquals(1, step1.getIndex());
        assertEquals(1, step1.getGeneration());
        assertEquals(FIRST, step1.getParent());

        List<Replicator> step2 = new ArrayList<Replicator>();

        step2.add(FIRST.replicate());
        step2.add(step1.replicate());

        assertEquals(2, step2.size());

        assertEquals(2, step2.get(0).getIndex());
        assertEquals(3, step2.get(1).getIndex());

        assertEquals(1, step2.get(0).getGeneration());
        assertEquals(2, step2.get(1).getGeneration());

        assertEquals(FIRST, step2.get(0).getParent());
        assertEquals(step1, step2.get(1).getParent());
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("jean.junit.ReplicatorTest");
    }
}
