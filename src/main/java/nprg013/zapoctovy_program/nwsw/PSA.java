package nprg013.zapoctovy_program.nwsw;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.ResourceBundle;

/**
 * The parent class to both of the algorithm classes, holding
 * fields and methods common to both.
 */
public abstract class PSA {
    /**
     * The first sequence.
     */
    String seq1;
    /**
     * The second sequence.
     */
    String seq2;
    SubstMatrix substMatrix;
    GapPenalty gapPenalty;
    /**
     * Collection of parents of each node in the table, used
     * for traceback over the filled in table.
     */
    Map<Integer,Integer> parents;
    /**
     * The table holding used to hold and calculate the values
     * of the algorithm computation.
     */
    double[][] table;
    /**
     * Results in form of the list in case there are multiple
     * equally correct alignments.
     * Due to the algorithm's character, stored sequences
     * are reversed.
     */
    List<StringBuilder[]> results;
    public PSA(String seq1path, String seq2path, String matrixPath, List<String> gapValue) throws IOException, SubstMatrix.MatrixDimensionException {
        LoadSeq seqs = new LoadSeq(seq1path, seq2path);
        this.seq1 = seqs.getFirst();
        this.seq2 = seqs.getSecond();
        this.substMatrix = new SubstMatrix(matrixPath);
        this.gapPenalty = new GapPenalty(gapValue);
    }
    public PSA(String seq1, String seq2, String matchValue, String mismatchValue, String gapValue) throws SubstMatrix.MatrixDimensionException {
        this.seq1 = seq1;
        this.seq2 = seq2;
        this.substMatrix = new SubstMatrix(matchValue, mismatchValue);
        this.gapPenalty = new GapPenalty(gapValue);
    }

    /**
     * Method for retrieving values of nodes represented as an int.
     *
     * @param position  Int representation of table position.
     * @return  Value held at the respective position in the table.
     */
    protected Double tableValue(int position) {

        return table[position / table[0].length][position % table[0].length];
    }

    /**
     * Method that determines the int value
     * of a position to the up and left from
     * the current one.
     *
     * @param position  The current position.
     * @return  The position adjacent in the up left direction.
     */
    protected int northWest(int position) {
        return position - 1 - table[0].length;
    }

    /**
     * Method that determines the int value
     * of a position to the left from the current one.
     *
     * @param position  The current position.
     * @return  The position to the left.
     */
    protected int left(int position) {
        return position - 1;
    }

    /**
     * Method that determines the int value
     * of a position directly above the current one.
     *
     * @param position  The current position.
     * @return  The position above the current one.
     */
    protected int up(int position) {
        return position - table[0].length;
    }

    /**
     * Calculates the position on the first sequence
     * corresponding to the provided position
     * in the table.
     *
     * @param current  The position in the table.
     * @return  Corresponding position of the sequence.
     */
    protected int seq1position(int current) {
        return (current / table[0].length) - 1;
    }

    /**
     * Calculates the position on the first sequence
     * corresponding to the provided position
     * in the table.
     *
     * @param current  The position in the table.
     * @return  Corresponding position of the sequence.
     */
    protected int seq2position(int current) {
        return (current % table[0].length) - 1;
    }

    /**
     * The required method that sets the table and determines
     * the calculation process according to the gap type.
     */
    abstract void calculate();

    /**
     * Method that visits each table position and assigns
     * the appropriate value using the method for affine
     * gap calculation.
     */
    void countAffValue() {
        int m = table.length;
        int n = table[0].length;
        for (int i = 0; i < m * n; i++) {
            if (i % n == 0) {continue;}
            if (i / n == 0) {continue;}
            assignAffValueAndParent(i - n - 1, i);
        }
        backtrack();
    }

    /**
     * Method that assigns the correct value
     * and determines paternity for each node
     * using affine gap penalties.
     *
     * @param i  The {@link #northWest(int position)} of current position.
     * @param i1  The current position.
     */
    abstract void assignAffValueAndParent(int i, int i1);

    /**
     * Method that visits each table position and assigns
     * the appropriate value using the method for linear
     * gap calculation.
     */
    void countLinValue() {
        int m = table.length;
        int n = table[0].length;
        //div n == radek
        //mod n == sloupec
        for (int i = 0; i < m * n; i++) {
            if (i % n == 0) {continue;}
            if (i / n == 0) {continue;}
            assignLinearValueAndParent(i - n - 1, i - 1, i - n, i);
        }
        backtrack();
    }

    /**
     * Method that assigns the correct value
     * and determines paternity for each node
     * using affine gap penalties.
     *
     * @param northWest  The {@link #northWest(int position)} of current position.
     * @param left  The node to the left of current position.
     * @param up  The node above the current position.
     * @param current  The current position.
     */
    abstract void assignLinearValueAndParent(int northWest, int left, int up, int current);

    /**
     * Method used to retrieve found alignments from the table.
     */
    abstract void backtrack();

    /**
     * Method used to print all the found results.
     */
    void printResults() {
        ResourceBundle bundle = ResourceBundle.getBundle("prompts");
        int counter = 1;
        for (StringBuilder[] result: results) {
            System.out.println(bundle.getString("option") + counter);
            System.out.println(result[0].reverse());
            System.out.println(result[1].reverse());
            counter++;
        }
    }

    /**
     * Method used to retrieve the results for further use.
     *
     * @return  List of pairs of sequences with the optimal alignment, which
     *          need to be calculated beforehand using the {@link #calculate()} method.
     */
    List<String[]> getResults() {
        List<String[]> output = new ArrayList<>();
        for (StringBuilder[] result: results) {
            output.add(new String[]{result[0].reverse().toString(), result[1].reverse().toString()});
        }
        return output;
    }
}
