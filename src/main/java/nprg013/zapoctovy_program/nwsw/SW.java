package nprg013.zapoctovy_program.nwsw;

import java.io.IOException;
import java.util.*;

/**
 * The class that implements the Smith-Waterman Pairwise Sequence Alignment algorithm.
 */
public class SW extends PSA {
    public SW(String seq1path, String seq2path, String matrixPath, List<String> gapValue) throws IOException, SubstMatrix.MatrixDimensionException {
        super(seq1path, seq2path, matrixPath, gapValue);
    }
    public SW(String seq1, String seq2, String matchValue, String mismatchValue, String gapValue) throws SubstMatrix.MatrixDimensionException {
        super(seq1, seq2, matchValue, mismatchValue, gapValue);
    }

    /**
     * Entry point for the calculation.
     * Sets the table and determines
     * next course according to the gap
     * penalty type.
     */
    @Override
    void calculate() {
        parents = new HashMap<>();
        table = new double[seq1.length() + 1][seq2.length() + 1];
        for (int i = 0; i < table.length; i++) {
            table[i][0] = 0.0;
        }
        Arrays.fill(table[0], 0.0);
        GapPenalty.Type gapType = gapPenalty.getType();
        switch (gapType) {
            case AFFINE:
                countAffValue();
                break;
            case LINEAR:
                countLinValue();
                break;
            case ERROR:
                throw new GapPenalty.GapPenaltyNumericalError("Gap penalty not properly initialized");
        }
    }

    /**
     * Method that determines the value of a node and potentially
     * assigns it a parent.
     * <p>
     * The value is determined as maximum of
     * 0, value of left node minus the gap penalty, value of the above
     * node minus the gap penalty, or value of the node to the up and
     * left incremented by the similarity score of corresponding symbols.
     * <p>
     * If the maximum turns out to be zero, no parents are assigned, otherwise
     * it is the node that was used to calculate this maximum.
     *
     * @param northWest  The {@link #northWest(int position)} of current position.
     * @param left  The node to the left of current position.
     * @param up  The node above the current position.
     * @param current  The current position.
     */
    @Override
    void assignLinearValueAndParent(int northWest, int left, int up, int current) {
        int parent = -1;
        Character seq1Char = seq1.charAt(seq1position(current));
        Character seq2Char = seq2.charAt(seq2position(current));
        Double comparisonValue = substMatrix.score(seq1Char, seq2Char);
        double addComparisonValue = tableValue(northWest) + comparisonValue;
        double endValue = Math.max(addComparisonValue, 0.0);
        double valueHolder = endValue;
        if (endValue > 0.0) {
            parent = northWest;
        }
        endValue = Math.max(endValue, tableValue(left) - gapPenalty.linearPart());
        if (endValue > valueHolder) {
            parent = left;
            valueHolder = endValue;
        }
        endValue = Math.max(endValue, tableValue(up) - gapPenalty.linearPart());
        if (endValue > valueHolder) {
            parent = up;
        }
        table[current / table[0].length][current % table[0].length] = endValue;
        if (parent >= 0) {
            parents.put(current, parent);
        }
    }

    /**
     * Method that determines the value of a node and potentially
     * assigns it a parent.
     * <p>
     * The value is determined as maximum of 0,
     * value of the node to the up and left
     * incremented by the similarity score of corresponding symbols,
     * maximal value  of some node above the current one minus the corresponding gap penalty,
     * or maximal value of some node to the left minus the corresponding gap penalty.
     * <p>
     * If the maximum turns out to be zero, no parents are assigned, otherwise
     * it is the node that was used to calculate this maximum.
     *
     * @param northWest  The node to the up left of the current one.
     * @param current The current node.
     */
    @Override
    void assignAffValueAndParent(int northWest, int current) {
        int parent = -1;
        int n = table[0].length;
        Character seq1Char = seq1.charAt(seq1position(current));
        Character seq2Char = seq2.charAt(seq2position(current));
        Double comparisonValue = substMatrix.score(seq1Char, seq2Char);
        double addComparisonValue = tableValue(northWest) + comparisonValue;
        double endValue = Math.max(addComparisonValue, 0.0);
        double valueHolder = endValue;
        if (endValue > 0.0) {
            parent = northWest;
        }
        double leftMaxValueCandidate = 0.0;
        int leftParentCandidate = -1;
        for (int i = 1; i < current % (n); i++) {
            if (tableValue(current - i) - gapPenalty.countAffine(i) > leftMaxValueCandidate) {
                leftMaxValueCandidate = tableValue(current-i) - gapPenalty.countAffine(i);
                leftParentCandidate = current - i;
            }
        }
        endValue = Math.max(endValue, leftMaxValueCandidate);
        if (endValue > valueHolder) {
            parent = leftParentCandidate;
            valueHolder = endValue;
        }
        double upMaxValueCandidate = 0.0;
        int upParentCandidate = -1;
        for (int i = 1; i < current / (n); i++) {
            if (tableValue(current - (n * i)) - gapPenalty.countAffine(i) > upMaxValueCandidate) {
                upMaxValueCandidate = tableValue(current - (n * i)) - gapPenalty.countAffine(i);
                upParentCandidate = current - (n * i);
            }
        }
        endValue = Math.max(endValue, upMaxValueCandidate);
        if (endValue > valueHolder) {
            parent = upParentCandidate;
        }
        table[current / n][current % n] = endValue;
        if (parent >= 0) {
            parents.put(current, parent);
        }
    }

    /**
     * Method used to retrieve the results from the completed table.
     * The process starts on the node with the highest value and follows
     * the succession of parents until it end on a node with a value of
     * zero, which is parentless.
     * <p>
     * Along this way, corresponding symbols are used to form the resulting
     * sequences from the end to the front.
     */
    @Override
    void backtrack() {
        int m = table.length;
        int n = table[0].length;
        results = new ArrayList<>();
        List<Integer> maxValuePosition = new ArrayList<>();
        Double currentMaxValue = 0.0;
        for (int i = 0; i < m * n; i++) {
            if (tableValue(i) > currentMaxValue) {
                currentMaxValue = tableValue(i);
            }
        }
        for (int i = 0; i < m * n; i++) {
            if (tableValue(i).equals(currentMaxValue)) {
                maxValuePosition.add(i);
            }
        }
        for (int current: maxValuePosition) {
            StringBuilder alignment1 = new StringBuilder();
            StringBuilder alignment2 = new StringBuilder();
            while(parents.containsKey(current)) {
                int parent = parents.get(current);
                if (parent == northWest(current)) {
                    alignment1.append(seq1.charAt(seq1position(current)));
                    alignment2.append(seq2.charAt(seq2position(current)));
                }
                else if (parent == left(current)) {
                    alignment1.append('_');
                    alignment2.append(seq2.charAt(seq2position(current)));
                }
                else if (parent == up(current)) {
                    alignment1.append(seq1.charAt(seq1position(current)));
                    alignment2.append('_');
                }
                else if (parent != left(current) && parent / n == current / n) {
                    while (current != parent) {
                        alignment1.append('_');
                        alignment2.append(seq2.charAt(seq2position(current)));
                        current--;
                    }
                }
                else if (parent != up(current) && parent % n == current % n) {
                    while (current != parent) {
                        alignment1.append(seq1.charAt(seq1position(current)));
                        alignment2.append('_');
                        current = current - n;
                    }
                }
                current = parent;
            }
            results.add(new StringBuilder[]{alignment1,alignment2});
        }
    }
}
