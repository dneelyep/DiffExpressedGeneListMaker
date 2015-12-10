import java.util.HashSet;
import java.util.Set;

/**
 * Represents a category into which a sample can be
 * placed for comparison with other samples.
 *
 * Each category has a set of valid values.
 *
 * An example category could be "Time", with values "0H", "24H', and "48H".
 */
public class SampleCategory {

    /**
     * A human-understandable name for the category.
     */
    public String Name;

    /// TODO Better docstring here.
    /**
     * A collection of all possible values that can belong to this category.
     */
    public HashSet<String> ValidValues;

    public SampleCategory(String name, HashSet<String> validValues) {
        this.Name = name;
        this.ValidValues = validValues;
    }

}
