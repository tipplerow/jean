
package jean.ensembl;

/**
 * Represents the unique Ensembl gene identifier.
 */
public final class EnsemblGeneID extends EnsemblID {
    private static final String LABEL_CODE = "gene:";

    private EnsemblGeneID(String key) {
        super(key);
    }

    /**
     * The canonical column name for gene identifiers in the header
     * line of data files to be analyzed by the {@code jean} library.
     */
    public static final String COLUMN_NAME = "Gene_ID";

    /**
     * Returns the Ensemble gene identifier for a given key string.
     *
     * @param key the key string.
     *
     * @return the Ensemble gene identifier for the given key string.
     */
    public static EnsemblGeneID instance(String key) {
        return new EnsemblGeneID(key);
    }

    /**
     * Extracts the gene key from an Ensembl record header line.
     *
     * @param headerLine the header line from an Ensembl record.
     *
     * @return the gene key contained in the given header line.
     *
     * @throws RuntimeException unless the header line contains a
     * properly formatted gene key.
     */
    public static EnsemblGeneID parseHeader(String headerLine) {
        return new EnsemblGeneID(parseHeader(headerLine, LABEL_CODE));
    }
}
