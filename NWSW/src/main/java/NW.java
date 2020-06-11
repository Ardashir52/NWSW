
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * The class that implements the Needleman-Wunsch Pairwise Sequence Alignment algorithm.
 */
public class NW extends PSA {
    Map<Integer,List<Integer>> parents;
    public NW(String seq1path, String seq2path, String matrixPath, List<String> gapValue) throws IOException, SubstMatrix.MatrixDimensionException {
        super(seq1path, seq2path, matrixPath, gapValue);
    }
    public NW(String seq1, String seq2, String match, String mismatch, String gapValue) throws SubstMatrix.MatrixDimensionException {
        super(seq1, seq2, match, mismatch, gapValue);
    }

    /**
     * Entry point for the calculation.
     * Sets the table and determines
     * next course according to the gap
     * penalty type. Also sets the parents for
     * the padding nodes as here everyone except
     * the first node must have a parent.
     */
    @Override
    void calculate() {
        parents = new HashMap<>();
        table = new Double[seq1.length() + 1][seq2.length() + 1];
        int n = table[0].length;
        String penaltyType = gapPenalty.getType();
        if (penaltyType.equals("affine")) {
            table[0][0] = 0.0;
            for (int i = 1; i < table.length; i++) {
                table[i][0] = (double) (-1) * gapPenalty.countAffine(i);
                parents.put((n * i), new ArrayList<>(List.of(n * (i - 1))));
            }
            for (int i = 1; i < n; i++) {
                table[0][i] = (double) (-1) * gapPenalty.countAffine(i);
                parents.put(i, new ArrayList<>(List.of(i - 1)));
            }
            countAffValue();
        }
        else if (penaltyType.equals("linear")) {
            for (int i = 0; i < table.length; i++) {
                table[i][0] = (double) i * (-1) * gapPenalty.linearPart();
                if (i != 0) {
                    parents.put(n * i, new ArrayList<>(List.of(n * (i - 1))));
                }
            }
            for (int i = 0; i < n; i++) {
                table[0][i] = (double) i * (-1) * gapPenalty.linearPart();
                if (i != 0) {
                    parents.put(i, new ArrayList<>(List.of(i - 1)));
                }
            }
            countLinValue();
        }
    }

    /**
     * This method fills in the table and assigns proper parents.
     * The value is calculated as a maximum of three values,
     * first the maximum of some node above the current one minus the gap penalty,
     * second the maximum of some node left of the current one minus the gap penalty,
     * third the value of the node to the up left plus the similarity score.
     * <p>
     * If any of these are equal, all of them are assign as a parent.
     *
     * @param northWest  The node to the up and left of the current position.
     * @param current  The current position.
     */
    @Override
    void assignAffValueAndParent(int northWest, int current) {
        List<Integer> parent = new ArrayList<>();
        int n = table[0].length;
        Character seq1Char = seq1.charAt(seq1position(current));
        Character seq2Char = seq2.charAt(seq2position(current));
        Double comparisonValue = substMatrix.score(seq1Char, seq2Char);
        double northWestValue = tableValue(northWest) + comparisonValue;
        double endValue = northWestValue;
        double leftValueCandidate = tableValue(current-1) - gapPenalty.countAffine(1);
        List<Integer> leftParentCandidate = new ArrayList<>();
        leftParentCandidate.add(current - 1);
        for (int i = 2; i < current % n; i++) {
            if (tableValue(current - i) - gapPenalty.countAffine(i) > leftValueCandidate) {
                leftValueCandidate = tableValue(current - i) - gapPenalty.countAffine(i);
                leftParentCandidate.clear();
                leftParentCandidate.add(current - i);
            }
            else if (tableValue(current - i) - gapPenalty.countAffine(i) == leftValueCandidate) {
                leftParentCandidate.add(current - i);
            }
        }
        double upValueCandidate = tableValue(current - n) - gapPenalty.countAffine(1);
        List<Integer> upParentCandidate = new ArrayList<>();
        upParentCandidate.add(current - n);
        for (int i = 2; i < current / n; i++) {
            if (tableValue(current - (n * i)) - gapPenalty.countAffine(i) > upValueCandidate) {
                upValueCandidate = tableValue(current - (n * i)) - gapPenalty.countAffine(i);
                upParentCandidate.clear();
                upParentCandidate.add(current- (n * i));
            }
            else if (tableValue(current - (n * i)) - gapPenalty.countAffine(i) == upValueCandidate) {
                upParentCandidate.add(current - (n * i));
            }
        }
        if (northWestValue >= upValueCandidate && northWestValue >= leftValueCandidate) {
            parent.add(northWest);
        }
        if (leftValueCandidate >= northWestValue && leftValueCandidate >= upValueCandidate) {
            parent.addAll(leftParentCandidate);
            endValue = leftValueCandidate;
        }
        if (upValueCandidate >= northWestValue && upValueCandidate >= leftValueCandidate) {
            parent.addAll(upParentCandidate);
            endValue = upValueCandidate;
        }
        table[current / n][current % n] = endValue;
        parents.put(current, parent);
    }

    /**
     * This method assigns each node its value and up to three parents.
     * The value is calculated as maximum of three values,
     * first the left node value minus the gap penalty,
     * second the up node value minus the gap penalty,
     * third the up left node plus the similarity score.
     * Each of these values equal to the maximum is kept as a parent.
     *
     * @param northWest  The {@link #northWest(int position)} of current position.
     * @param left  The node to the left of current position.
     * @param up  The node above the current position.
     * @param current  The current position.
     */
    @Override
    void assignLinearValueAndParent(int northWest, int left, int up, int current) {
        List<Integer> parent = new ArrayList<>();
        Character seq1Char = seq1.charAt(seq1position(current));
        Character seq2Char = seq2.charAt(seq2position(current));
        Double comparisonValue = substMatrix.score(seq1Char, seq2Char);
        double northWestValue = tableValue(northWest) + comparisonValue;
        double leftValue = tableValue(left) - gapPenalty.linearPart();
        double upValue = tableValue(up) - gapPenalty.linearPart();
        double endValue = northWestValue;
        if (northWestValue >= leftValue && northWestValue >= upValue) {
            parent.add(northWest);
        }
        if (leftValue >= northWestValue && leftValue >= upValue) {
            parent.add(left);
            endValue = leftValue;
        }
        if (upValue >= northWestValue && upValue >= leftValue) {
            parent.add(up);
            endValue = upValue;
        }
        table[current / table[0].length][current % table[0].length] = endValue;
        parents.put(current, parent);
    }

    /**
     * The method for retrieving the results.
     * The process starts on the last node with empty sequences
     * and is immediately passed to the {@link #backtrackTree(StringBuilder, StringBuilder, int)} method.
     */
    @Override
    void backtrack() {
        results = new ArrayList<>();
        int m = table.length;
        int n = table[0].length;
        int start = m * n -1;
        backtrackTree(new StringBuilder(), new StringBuilder(), start);
    }

    /**
     * A method that finds all the possible results
     * by exploring each parent option of each node.
     *
     * @param sequence1  The first sequence part so far built.
     * @param sequence2  The second sequence part so far built.
     * @param position  The current position.
     */
    void backtrackTree(StringBuilder sequence1, StringBuilder sequence2, int position) {
        int n = table[0].length;
        if (position == 0) {
            results.add(new StringBuilder[]{new StringBuilder(sequence1), new StringBuilder(sequence2)});
        }
        else {
            List<Integer> listOfParents = parents.get(position);
            for (Integer parent: listOfParents) {
                if (parent == northWest(position)) {
                    StringBuilder newSeq1 = new StringBuilder(sequence1);
                    StringBuilder newSeq2 = new StringBuilder(sequence2);
                    newSeq1.append(seq1.charAt(seq1position(position)));
                    newSeq2.append(seq2.charAt(seq2position(position)));
                    backtrackTree(newSeq1, newSeq2, parent);
                }
                else if (parent / n == position / n) {
                    StringBuilder newSeq1 = new StringBuilder(sequence1);
                    StringBuilder newSeq2 = new StringBuilder(sequence2);
                    while (position != parent) {
                        newSeq1.append('_');
                        newSeq2.append(seq2.charAt(seq2position(position)));
                        position--;
                    }
                    backtrackTree(newSeq1, newSeq2, parent);
                }
                else if (parent % n == position % n) {
                    StringBuilder newSeq1 = new StringBuilder(sequence1);
                    StringBuilder newSeq2 = new StringBuilder(sequence2);
                    while (position != parent) {
                        newSeq1.append(seq1.charAt(seq1position(position)));
                        newSeq2.append('_');
                        position = position - n;
                    }
                    backtrackTree(newSeq1, newSeq2, parent);
                }
            }
        }
    }
}
