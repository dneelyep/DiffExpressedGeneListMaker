import java.io.File;

/**
 * Represents a set of parameters that should be applied
 * when generating a set of differentially-expressed genes.
 */
public class GeneListParameters {

    /**
     * The directory to retrieve CEL files from.
     */
    public File CELFileDirectory;

    public GeneListParameters(File celFileDirectory) {
        this.CELFileDirectory = celFileDirectory;
    }
}