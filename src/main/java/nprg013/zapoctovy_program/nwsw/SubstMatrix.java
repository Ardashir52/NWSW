package nprg013.zapoctovy_program.nwsw;

import java.io.*;
import java.util.*;

/**
 * This class is used to load and keep the values of substitution matrices,
 * that is the values of matches and mismatches between the symbols in the
 * sequences considered.
 */
public class SubstMatrix {
    /**
     * List of all symbols ordered by their position in the source matrix.
     */
    private List<Character> keys;
    /**
     * The numerical values of the source matrix represented by a 2D table
     * of doubles.
     */
    private final Double[][] values;
    /**
     * Determines, whether the matrix was constructed under the interactive mode.
     */
    private final boolean interactive;

    /**
     * Reads a matrix from an input file and parses it to the respective
     * fields.
     * <p>
     * The 0,0 position of the original matrix will be ignored, as it is
     * irrelevant and indeed customarily left blank in common usage.
     *
     * @param path  The path to file with the scoring matrix.
     * @throws IOException  This is thrown if the path provided is not valid.
     * @throws MatrixDimensionException  This is thrown in case the matrix is not symmetrical
     *                                   or the values cannot be parsed.
     */
    public SubstMatrix(String path) throws IOException, MatrixDimensionException {
        this.interactive = false;
        File f = new File(path);
        List<String[]> vstupniTabulka = new ArrayList<>();
        try (BufferedReader bReader = new BufferedReader(new FileReader(f));) {
            String line;
            while ((line = bReader.readLine()) != null) {
                String[] elements = line.trim().split("\\s+");
                vstupniTabulka.add(elements);
            }
            keys = new ArrayList<>();
            String[] keysArr;
            try {
                if (vstupniTabulka.get(0).length == vstupniTabulka.get(1).length) {
                    keysArr = Arrays.copyOfRange(vstupniTabulka.get(0), 1, vstupniTabulka.get(0).length);
                } else {
                    keysArr = Arrays.copyOfRange(vstupniTabulka.get(0), 0, vstupniTabulka.get(0).length);
                }
                for (String key : keysArr) {
                    if (key.length() == 1) {
                        keys.add(key.charAt(0));
                        continue;
                    }
                    keys.add(oneCharCode(key));
                }
                values = new Double[keys.size()][keys.size()];
                for (int i = 1; i < vstupniTabulka.size(); i++) {
                    for (int j = 1; j < vstupniTabulka.get(i).length; j++) {
                        values[i - 1][j - 1] = Double.parseDouble(vstupniTabulka.get(i)[j]);
                    }
                }
            } catch (Exception e) {
                throw new MatrixDimensionException("Wrong matrix format!");
            }
        }
    }

    /**
     * The matrix constructor adapted for the interactive mode.
     *
     * @param match  The value used for matches
     * @param mismatch  The value used if the two characters under consideration are different.
     * @throws MatrixDimensionException  This is thrown in case the parameters cannot be parsed
     *                                   to double.
     */
    public SubstMatrix(String match, String mismatch) throws MatrixDimensionException {
        this.interactive = true;
        try {
            values = new Double[1][2];
            values[0][0] = Double.parseDouble(match);
            values[0][1] = Double.parseDouble(mismatch);
        }
        catch (Exception e) {
            throw new MatrixDimensionException("Wrong match/mismatch values");
        }
    }

    /**
     * This method returns the similarity score of two symbols.
     * In the case of interactive mode just determines whether
     * they are identical or not and assigns score accordingly.
     *
     * @param a  First character for comparison
     * @param b  Second character for comparison
     * @return  The table value of similarity of the two parameters.
     */
    public Double score(Character a, Character b) {
        if (interactive) {
            if (a.equals(b)) {
                return values[0][0];
            }
            else {
                return values[0][1];
            }
        }
        else {
            int keyA = keys.indexOf(a);
            int keyB = keys.indexOf(b);
            return values[keyA][keyB];
        }
    }

    /**
     * Conversion of the three letter codes of amino acids to corresponding
     * one letter codes storable in the {@link #keys} field.
     *
     * @param threeCharCode  The three letter amino acid code.
     * @return  The corresponding one letter amino acid code.
     */
    private Character oneCharCode(String threeCharCode) {
        switch (threeCharCode.toLowerCase()) {
            case "ala":
                return 'A';
            case "arg":
                return 'R';
            case "asn":
                return 'N';
            case "asp":
                return 'D';
            case "asx":
                return 'B';
            case "cys":
                return 'C';
            case "glu":
                return 'E';
            case "gln":
                return 'Q';
            case "glx":
                return 'Z';
            case "gly":
                return 'G';
            case "his":
                return 'H';
            case "ile":
                return 'I';
            case "leu":
                return 'L';
            case "lys":
                return 'K';
            case "met":
                return 'M';
            case "phe":
                return 'F';
            case "pro":
                return 'P';
            case "ser":
                return 'S';
            case "thr":
                return 'T';
            case "trp":
                return 'W';
            case "tyr":
                return 'Y';
            case "val":
                return 'V';
            default:
                return '?';
        }
    }

    /**
     * This exception is raised in the case that the matrix is not symmetrical mainly,
     * but would also apply if the input values cannot be parsed.
     */
    public static class MatrixDimensionException extends Exception {
        public MatrixDimensionException(String s) {
            super(s);
        }
    }
}
