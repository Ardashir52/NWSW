package nprg013.zapoctovy_program.nwsw;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * This rather simple class has the loading of the examined sequences for its task.
 */
public class LoadSeq {
    /**
     * One of the sequences to compare and sequence info extracted from
     * the comments lines of the FASTA format.
     */
    private final SequenceData seq1;
    /**
     * One of the sequences to compare and sequence info extracted from
     * the comments lines of the FASTA format.
     */
    private final SequenceData seq2;

    /**
     * Loads the two sequences and potential comments.
     *
     * @param path1  The path to the file with the first sequence.
     * @param path2  The path to the file with the second sequence.
     * @throws IOException  Thrown in case at least one of the paths is not valid.
     */
    public LoadSeq(String path1, String path2) throws IOException {
        File f = new File(path1);
        try (BufferedReader reader = new BufferedReader(new FileReader(f))) {
            seq1 = readSequence(reader);
        }
        f = new File(path2);
        try (BufferedReader reader = new BufferedReader(new FileReader(f))) {
            seq2 = readSequence(reader);
        }
    }

    public LoadSeq(InputStream in1, InputStream in2) throws IOException {
        try (BufferedReader firstReader = new BufferedReader(new InputStreamReader(in1));
             BufferedReader secondReader = new BufferedReader(new InputStreamReader(in2))) {
            seq1 = readSequence(firstReader);
            seq2 = readSequence(secondReader);
        }
    }

    /**
     * Returns the first sequence.
     *
     * @return  The first sequence.
     */
    public String getFirst() { return seq1.getValue(); }

    /**
     * Returns the second sequence.
     *
     * @return  The second sequence.
     */
    public String getSecond() {
        return seq2.getValue();
    }

    /**
     * Returns facultative info about the first sequence.
     *
     * @return  Possibly empty list of strings representing the lines of commentary.
     */
    public List<String> firstInfo() {
        return seq1.getInfo();
    }

    /**
     * Returns facultative info about the first sequence.
     *
     * @return  Possibly empty list of strings representing the lines of commentary.
     */
    public List<String> secondInfo() {
        return seq2.getInfo();
    }

    /**
     * Method reads sequence and its info from a file.
     *
     * @param bReader  Initialized reader object.
     * @return  {@link SequenceData} type object.
     * @throws IOException  In case of file related error.
     */
    SequenceData readSequence(BufferedReader bReader) throws IOException {
        StringBuilder seqBuilder = new StringBuilder();
        List<String> info = new ArrayList<>();
        String line;
        while ((line = bReader.readLine()) != null) {
            if (line.trim().charAt(0) == '>') {
                info.add(line);
            }
            else {
                seqBuilder.append(line.trim());
            }
        }
        return new SequenceData(seqBuilder.toString(), info);
    }

    /**
     * Class designated to store the sequence data.
     */
    static class SequenceData {
        String sequenceValue;
        List<String> sequenceInfo;
        public SequenceData(String value, List<String> info) {
            this.sequenceInfo = info;
            this.sequenceValue = value;
        }

        /**
         * Returns the value of the sequence.
         *
         * @return String value of sequence.
         */
        String getValue() {
            return this.sequenceValue;
        }

        /**
         * Returns the comment part of the sequence.
         *
         * @return  Comment part of sequence as a List of String
         */
        List<String> getInfo() {
            return this.sequenceInfo;
        }
    }
}
