
package jean.ensembl;

import java.io.File;
import java.io.FileFilter;

import jam.app.JamEnv;
import jam.app.JamProperties;
import jam.lang.JamException;

/**
 * Resolves path names to Ensembl data files.
 */
public final class EnsemblLocator {
    /**
     * Environment variable that defines the absolute path name for
     * the directory that contains Ensembl chromosome FASTA files.
     * If the system property {@code jean.ensembl.genomeDir} is also
     * defined, it will override the environment variable.
     */
    public static final String GENOME_DIR_ENV = "JEAN_ENSEMBL_GENOME_DIR";

    /**
     * System property that defines the absolute path name for the
     * directory that contains Ensembl chromosome FASTA files.  If
     * not defined, environment variable {@code GENOME_DIR_ENV} will
     * be used by default.
     */
    public static final String GENOME_DIR_PROPERTY = "jean.ensembl.genomeDir";

    /**
     * Returns the file that contains the nucleotide sequence for a
     * specific human chromosome.
     *
     * @param code the code (number or letter) of the chromosome to
     * retrieve.
     *
     * @return the file that contains the nucleotide sequence for the
     * specified human chromosome.
     *
     * @throws RuntimeException unless a FASTA file for the specified
     * chromosome exists in the genome directory.
     */
    public static File resolveChromosomeFile(String code) {
        File[] files = resolveGenomeDir().listFiles(chromosomeFileFilter(code));

        if (files.length == 1)
            return files[0];
        else
            throw JamException.runtime("Could not find nucleotide sequence for chromosome [%s].", code);
    }

    private static FileFilter chromosomeFileFilter(String code) {
        return new FileFilter() {
            @Override public boolean accept(File file) {
                String baseName = file.getName();
                String speciesFragment = "Homo_sapiens.GRCh";
                String chromosomeFragment = ".dna.chromosome." + code + ".fa";

                return baseName.startsWith(speciesFragment)
                    && baseName.contains(chromosomeFragment);
            }
        };
    }

    /**
     * Returns the directory that contains the Ensembl chromosome
     * FASTA files.
     *
     * @return the directory that contains the Ensembl chromosome
     * FASTA files.
     */
    public static File resolveGenomeDir() {
        return new File(resolveGenomeDirName());
    }

    /**
     * Returns the name of the directory that contains the Ensembl
     * chromosome FASTA files.
     *
     * @return the name of the directory that contains the Ensembl
     * chromosome FASTA files.
     */
    public static String resolveGenomeDirName() {
        if (JamProperties.isSet(GENOME_DIR_PROPERTY))
            return JamProperties.getRequired(GENOME_DIR_PROPERTY);
        else
            return JamEnv.getRequired(GENOME_DIR_ENV);
    }
}
