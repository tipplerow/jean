
package jean.ensembl;

/**
 * Represents the unique Ensembl transcript identifier.
 */
public final class EnsemblTranscriptID extends EnsemblID {
    private static final String LABEL_CODE = "transcript:";

    private EnsemblTranscriptID(String key) {
        super(key);
    }

    /**
     * The canonical column name for transcript identifiers in the
     * header line of data files to be analyzed by the {@code jean}
     * library.
     */
    public static final String COLUMN_NAME = "Transcript_ID";

    /**
     * Returns the Ensemble transcript identifier for a given key
     * string.
     *
     * @param key the key string.
     *
     * @return the Ensemble transcript identifier for the given key
     * string.
     */
    public static EnsemblTranscriptID instance(String key) {
        return new EnsemblTranscriptID(key);
    }

    /**
     * Extracts the transcript key from an Ensembl record header line.
     *
     * @param headerLine the header line from an Ensembl record.
     *
     * @return the transcript key contained in the given header line.
     *
     * @throws RuntimeException unless the header line contains a
     * properly formatted transcript key.
     */
    public static EnsemblTranscriptID parseHeader(String headerLine) {
        return new EnsemblTranscriptID(parseHeader(headerLine, LABEL_CODE));
    }
}
