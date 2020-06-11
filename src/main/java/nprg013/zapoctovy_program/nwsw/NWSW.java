package nprg013.zapoctovy_program.nwsw;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.ResourceBundle;

/**
 * A class that performs Pairwise Sequence Alignment by the means of two algorithms,
 * namely Needleman-Wunsch and Smith-Waterman.
 */
public class NWSW {
    /**
     * The entry point of the program. Absence of arguments starts the interactive mode, whereas
     * expected argument count for basic function is four or five, depending on the type of gap penalty.
     * All other argument counts will result in help being displayed.
     *
     * @param argv  Array of arguments.
     * @throws IOException  Thrown on failure to read input.
     */
    public static void main(String[] argv) throws IOException {
        if (argv.length > 3 && argv.length < 6) {
            try {
                List<String> gapArgs = new ArrayList<>();
                gapArgs.add(argv[3]);
                if (argv.length == 5) {
                    gapArgs.add(argv[4]);
                }
                boolean swRatherThanNW = userChoice();
                if (swRatherThanNW) {
                    SW swRun = new SW(argv[0], argv[1], argv[2], gapArgs);
                    swRun.calculate();
                    swRun.printResults();
                }
                else {
                    NW nwRun = new NW(argv[0], argv[1], argv[2], gapArgs);
                    nwRun.calculate();
                    nwRun.printResults();
                }
            } catch (GapPenalty.GapPenaltyNumericalError error) {
                error.printStackTrace();
                help();
            } catch (Exception error) {
                error.printStackTrace();
            }
        }
        else if (argv.length != 0){
            help();
        }
        else {
            interactiveRun();
        }

    }

    /**
     * Describes intended argument usage.
     */
    static void help() {
        ResourceBundle bundle = ResourceBundle.getBundle("prompts");
        System.out.println(bundle.getString("help1"));
        System.out.println(bundle.getString("help2"));
        System.out.println(bundle.getString("help3"));
    }

    /**
     * Asks the user about their preference of PSA algorithm.
     *
     * @return  True if user chooses SW, false for NW.
     */
    static boolean userChoice() throws IOException {
        ResourceBundle bundle = ResourceBundle.getBundle("prompts");
        System.out.println(bundle.getString("choice1"));
        System.out.println(bundle.getString("choice2"));
        System.out.println(bundle.getString("choice3"));
        BufferedReader bReader = new BufferedReader(new InputStreamReader(System.in));
        String userAnswer = bReader.readLine();
        while (!userAnswer.toUpperCase().equals("A") && !userAnswer.toUpperCase().equals("B")) {
            System.out.println(bundle.getString("choice4"));
            userAnswer = bReader.readLine();
        }
        return userAnswer.toUpperCase().equals("B");
    }

    /**
     * Method that realizes the interactive mode.
     *
     * @throws IOException On input reading error.
     */
    static void interactiveRun() throws IOException {
        ResourceBundle bundle = ResourceBundle.getBundle("prompts");
        System.out.println(bundle.getString("inter1"));
        BufferedReader bReader = new BufferedReader(new InputStreamReader(System.in));
        String vstup = bReader.readLine();
        while (!vstup.toLowerCase().equals("help") && !vstup.toLowerCase().equals("inter") ) {
            System.out.println(bundle.getString("inter2"));
            vstup = bReader.readLine();
        }
        switch (vstup) {
            case "help":
                help();
                break;
            case "inter":
                System.out.println(bundle.getString("inter3"));
                String seq1 = bReader.readLine();
                System.out.println(bundle.getString("inter4"));
                String seq2 = bReader.readLine();
                System.out.println(bundle.getString("inter5"));
                String match = bReader.readLine();
                System.out.println(bundle.getString("inter6"));
                String mismatch = bReader.readLine();
                System.out.println(bundle.getString("inter7"));
                String gp = bReader.readLine();
                System.out.println(bundle.getString("inter8"));
                System.out.println(bundle.getString("inter9"));
                String option = bReader.readLine();
                while (!option.toLowerCase().equals("s") && !option.toLowerCase().equals("n")) {
                    System.out.println(bundle.getString("inter10"));
                    option = bReader.readLine();
                }
                if (option.toLowerCase().equals("s")) {
                    try {
                        SW interSW = new SW(seq1, seq2, match, mismatch, gp);
                        interSW.calculate();
                        interSW.printResults();
                    } catch (SubstMatrix.MatrixDimensionException | GapPenalty.GapPenaltyNumericalError error) {
                        error.printStackTrace();
                    }
                } else {
                    try {
                        NW interNW = new NW(seq1, seq2, match, mismatch, gp);
                        interNW.calculate();
                        interNW.printResults();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
        }
    }

    /**
     * Compares two sequences by the means of the Needleman-Wunsch algorithm using affine gap penalty.
     *
     * @param seq1path  Path to file with first sequence.
     * @param seq2path  Path to file with second sequence.
     * @param mxPath  Path to file with the scoring matrix.
     * @param gpLinear  Value of the linear part of the gap penalty.
     * @param gpAff  Value of the affine part of the gap penalty.
     * @return  A list of all the optimal alignments.
     * @throws IOException  If some of the provided paths is not valid.
     * @throws SubstMatrix.MatrixDimensionException  If the provided matrix is not valid.
     */
    public List<String[]> needlemanWunsch(String seq1path, String seq2path, String mxPath, double gpLinear, double gpAff) throws IOException, SubstMatrix.MatrixDimensionException {
        NW nw = new NW(seq1path, seq2path, mxPath, new ArrayList<>(List.of(Double.toString(gpLinear), Double.toString(gpAff))));
        nw.calculate();
        return nw.getResults();
    }

    /**
     * Compares two sequences by the means of the Needleman-Wunsch algorithm using linear gap penalty.
     *
     * @param seq1path  Path to file with first sequence.
     * @param seq2path  Path to file with second sequence.
     * @param mxPath  Path to file with the scoring matrix.
     * @param gpLinear  Value of the linear part of the gap penalty.
     * @return  A list of all the optimal alignments.
     * @throws IOException  If some of the provided paths is not valid.
     * @throws SubstMatrix.MatrixDimensionException  If the provided matrix is not valid.
     */
    public List<String[]> needlemanWunsch(String seq1path, String seq2path, String mxPath, double gpLinear) throws IOException, SubstMatrix.MatrixDimensionException {
        NW nw = new NW(seq1path, seq2path, mxPath, new ArrayList<>(List.of(Double.toString(gpLinear))));
        nw.calculate();
        return nw.getResults();
    }

    /**
     * Compares two sequences by the means of the Smith-Waterman algorithm using affine gap penalty.
     *
     * @param seq1path  Path to file with first sequence.
     * @param seq2path  Path to file with second sequence.
     * @param mxPath  Path to file with the scoring matrix.
     * @param gpLinear  Value of the linear part of the gap penalty.
     * @param gpAff  Value of the affine part of the gap penalty.
     * @return  A list of all the optimal alignments.
     * @throws IOException  If some of the provided paths is not valid.
     * @throws SubstMatrix.MatrixDimensionException  If the provided matrix is not valid.
     */
    public List<String[]> smithWaterman(String seq1path, String seq2path, String mxPath, double gpLinear, double gpAff) throws IOException, SubstMatrix.MatrixDimensionException {
        SW sw = new SW(seq1path, seq2path, mxPath, new ArrayList<>(List.of(Double.toString(gpLinear), Double.toString(gpAff))));
        sw.calculate();
        return sw.getResults();
    }

    /**
     * Compares two sequences by the means of the Smith-Waterman algorithm using linear gap penalty.
     *
     * @param seq1path  Path to file with first sequence.
     * @param seq2path  Path to file with second sequence.
     * @param mxPath  Path to file with the scoring matrix.
     * @param gpLinear  Value of the linear part of the gap penalty.
     * @return  A list of all the optimal alignments.
     * @throws IOException  If some of the provided paths is not valid.
     * @throws SubstMatrix.MatrixDimensionException  If the provided matrix is not valid.
     */
    public List<String[]> smithWaterman(String seq1path, String seq2path, String mxPath, double gpLinear) throws IOException, SubstMatrix.MatrixDimensionException {
        SW sw = new SW(seq1path, seq2path, mxPath, new ArrayList<>(List.of(Double.toString(gpLinear))));
        sw.calculate();
        return sw.getResults();
    }
}
