
package jean.ensembl;

import jam.util.RegexUtil;

/**
 * Represents the unique Ensembl protein identifier.
 */
public final class EnsemblProteinID extends EnsemblID {
    private EnsemblProteinID(String key) {
        super(key);
    }

    /**
     * Returns the Ensemble protein identifier for a given key string.
     *
     * @param key the key string.
     *
     * @return the Ensemble protein identifier for the given key string.
     */
    public static EnsemblProteinID instance(String key) {
        return new EnsemblProteinID(key);
    }

    /**
     * Extracts the protein key from an Ensembl record header line.
     *
     * @param headerLine the header line from an Ensembl record.
     *
     * @return the protein key contained in the given header line.
     *
     * @throws RuntimeException unless the header line contains a
     * properly formatted protein key.
     */
    public static EnsemblProteinID parseHeader(String headerLine) {
        return parseKey(RegexUtil.split(HEADER_FIELD_DELIM, headerLine)[0]);
    }

    /**
     * Extracts the protein key from a FASTA record key (by stripping
     * the version number).
     *
     * @param fastaKey the FASTA key from an Ensembl record.
     *
     * @return the protein key contained in the given FASTA key.
     */
    public static EnsemblProteinID parseKey(String fastaKey) {
        return instance(stripVersion(fastaKey));
    }
}
