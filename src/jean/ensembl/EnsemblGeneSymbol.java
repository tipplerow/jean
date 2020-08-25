
package jean.ensembl;

import jean.hugo.HugoSymbol;

/**
 * Represents the Ensembl gene symbol (HUGO) identifier.
 */
public final class EnsemblGeneSymbol extends EnsemblID {
    private static final String LABEL_CODE = "gene_symbol:";

    private EnsemblGeneSymbol(String key) {
        super(key);
    }

    /**
     * Extracts the HUGO symbol from an Ensembl record header line.
     *
     * @param headerLine the header line from an Ensembl record.
     *
     * @return the HUGO symbol contained in the given header line.
     *
     * @throws RuntimeException unless the header line contains a
     * properly formatted HUGO symbol.
     */
    public static HugoSymbol parseHeader(String headerLine) {
        return HugoSymbol.instance(parseHeader(headerLine, LABEL_CODE));
    }
}
