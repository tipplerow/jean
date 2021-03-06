
package jean.neo;

import java.util.Comparator;

import jam.io.Delimiter;
import jam.lang.JamException;
import jam.math.IntRange;
import jam.report.LineBuilder;

import jean.hugo.HugoSymbol;
import jean.peptide.Peptide;
import jean.tcga.TumorBarcode;
import jean.tcga.TumorGeneRecord;

/**
 * Associates a neo-peptide with the tumor sample, gene, and
 * self-peptide from which it originated.
 */
public final class PeptidePairRecord extends TumorGeneRecord {
    private final IntRange peptideRange;
    private final PeptidePair peptidePair;

    private PeptidePairRecord(TumorBarcode tumorBarcode,
                              HugoSymbol   hugoSymbol,
                              IntRange     peptideRange,
                              PeptidePair  peptidePair) {
        super(tumorBarcode, hugoSymbol);
        this.peptideRange = peptideRange;
        this.peptidePair  = peptidePair;
    }

    /**
     * The standard delimiter for flat files containing peptide pair
     * records.
     */
    public static final Delimiter DELIM = Delimiter.TAB;

    /**
     * A comparator that orders records by tumor barcode first, HUGO
     * symbol second, and peptide range third.
     */
    public static final Comparator<PeptidePairRecord> COMPARATOR =
        new Comparator<PeptidePairRecord>() {
            @Override public int compare(PeptidePairRecord rec1, PeptidePairRecord rec2) {
                int barcodeSymbolCmp = BARCODE_SYMBOL_COMPARATOR.compare(rec1, rec2);

                if (barcodeSymbolCmp != 0)
                    return barcodeSymbolCmp;
                else
                    return IntRange.BOUND_COMPARATOR.compare(rec1.peptideRange, rec2.peptideRange);
            }
        };

    /**
     * Returns a peptide pair record with fixed components.
     *
     * @param tumorBarcode the tumor where the mutation occurred.
     *
     * @param hugoSymbol the HUGO symbol of the mutated gene.
     *
     * @param peptideRange the unit-offset range of the amino acid
     * positions in the peptide fragments.
     *
     * @param peptidePair the peptide pair resulting from the
     * mutation.
     *
     * @return the peptide pair record with the specified
     * components.
     */
    public static PeptidePairRecord instance(TumorBarcode tumorBarcode,
                                             HugoSymbol   hugoSymbol,
                                             IntRange     peptideRange,
                                             PeptidePair  peptidePair) {
        return new PeptidePairRecord(tumorBarcode, hugoSymbol, peptideRange, peptidePair);
    }

    /**
     * Returns a peptide pair record with fixed components.
     *
     * @param tumorBarcode the tumor where the mutation occurred.
     *
     * @param hugoSymbol the HUGO symbol of the mutated gene.
     *
     * @param peptideRange the unit-offset range of the amino acid
     * positions in the peptide fragments.
     *
     * @param selfPeptide the germline self-peptide.
     *
     * @param neoPeptide the neo-peptide generated by somatic
     * mutation.
     *
     * @return the peptide pair record with the specified
     * components.
     */
    public static PeptidePairRecord instance(TumorBarcode tumorBarcode,
                                             HugoSymbol   hugoSymbol,
                                             IntRange     peptideRange,
                                             SelfPeptide  selfPeptide,
                                             NeoPeptide   neoPeptide) {
        return instance(tumorBarcode, hugoSymbol, peptideRange,
                        PeptidePair.instance(selfPeptide, neoPeptide));
    }                                            

    /**
     * Returns the header line for flat files containing peptide pair
     * records.
     *
     * @return the header line for flat files containing peptide pair
     * records.
     */
    public static String header() {
        LineBuilder builder = new LineBuilder(DELIM);

        builder.append(TumorBarcode.COLUMN_NAME);
        builder.append(HugoSymbol.COLUMN_NAME);
        builder.append("Range_Lower");
        builder.append("Range_Upper");
        builder.append("Self_Peptide");
        builder.append("Neo_Peptide");

        return builder.toString();
    }

    /**
     * Creates a new peptide pair record by parsing a delimited line
     * from a flat file.
     *
     * @param line the line to parse.
     *
     * @return the peptide pair record encoded in the specified line.
     *
     * @throws RuntimeException unless the line contains a properly
     * formatted peptide pair record.
     */
    public static PeptidePairRecord parse(String line) {
        String[] fields = DELIM.split(line, 6);

        TumorBarcode tumorBarcode = TumorBarcode.instance(fields[0]);
        HugoSymbol   hugoSymbol   = HugoSymbol.instance(fields[1]);
        int          rangeLower   = Integer.parseInt(fields[2]);
        int          rangeUpper   = Integer.parseInt(fields[3]);
        SelfPeptide  selfPeptide  = SelfPeptide.instance(fields[4]);
        NeoPeptide   neoPeptide   = NeoPeptide.instance(fields[5]);

        return instance(tumorBarcode, hugoSymbol,
                        IntRange.instance(rangeLower, rangeUpper),
                        PeptidePair.instance(selfPeptide, neoPeptide));
    }

    /**
     * Formats this record for output to a delimited flat file.
     *
     * @return a string containing the formatted text.
     */
    public String format() {
        LineBuilder builder = new LineBuilder(DELIM);

        builder.append(tumorBarcode.getKey());
        builder.append(hugoSymbol.getKey());
        builder.append(peptideRange.lower());
        builder.append(peptideRange.upper());
        builder.append(peptidePair.self().formatString());
        builder.append(peptidePair.neo().formatString());

        return builder.toString();
    }

    /**
     * Returns the unit-offset range of the amino acid positions in
     * the peptide fragments.
     *
     * @return the unit-offset range of the amino acid positions in
     * the peptide fragments.
     */
    public IntRange getPeptideRange() {
        return peptideRange;
    }

    /**
     * Returns the neo-peptide derived from the self-peptide by
     * somatic mutation.
     *
     * @return the neo-peptide derived from the self-peptide by
     * somatic mutation.
     */
    public Peptide getNeoPeptide() {
        return peptidePair.neo();
    }

    /**
     * Returns the self-peptide derived from the germline genome.
     *
     * @return the self-peptide derived from the germline genome.
     */
    public Peptide getSelfPeptide() {
        return peptidePair.self();
    }

    @Override public String toString() {
        return String.format("PeptidePairRecord(%s, %s, [%d, %d]: %s => %s)",
                             tumorBarcode.getKey(),
                             hugoSymbol.getKey(),
                             peptideRange.lower(),
                             peptideRange.upper(),
                             getSelfPeptide().formatString(),
                             getNeoPeptide().formatString());
    }
}
