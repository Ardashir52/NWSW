import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.Objects;

public class NWSWTest {
    @Test
    public void loadSeqTest() throws IOException {
        InputStream in1 = NWSWTest.class.getClassLoader().getResourceAsStream("sequences/P59594.fasta.txt");
        InputStream in2 = NWSWTest.class.getClassLoader().getResourceAsStream("sequences/P59594.fasta.txt");
        LoadSeq sequence = new LoadSeq(in1, in2);
        Assert.assertEquals(1255, sequence.getFirst().length());
        Assert.assertTrue(sequence.secondInfo().get(0).contains("coronavirus"));
    }
    @Test
    public void loadMXTest() throws IOException, SubstMatrix.MatrixDimensionException {
        String mxPath = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("matrices/PAM250.txt")).getPath();
        SubstMatrix matice = new SubstMatrix(mxPath);
        Assert.assertEquals(3.5, matice.score('V', 'V'), 0.0);
    }
    @Test
    public void swTestsSimple1() throws IOException, SubstMatrix.MatrixDimensionException {
        String seq1path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/test1a.txt")).getPath();
        String seq2path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/test1b.txt")).getPath();
        String mxPath = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("matrices/testMX.txt")).getPath();
        NWSW swTest = new NWSW();
        List<String[]> results = swTest.smithWaterman(seq1path, seq2path, mxPath, 2);
        Assert.assertEquals(1, results.size());
    }
    @Test
    public void swTestSimple2() throws IOException, SubstMatrix.MatrixDimensionException {
        String seq1path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/test1a.txt")).getPath();
        String seq2path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/test1b.txt")).getPath();
        String mxPath = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("matrices/testMX.txt")).getPath();
        NWSW swTest = new NWSW();
        List<String[]> results = swTest.smithWaterman(seq1path, seq2path, mxPath, 2);
        String first = results.get(0)[0];
        String second = results.get(0)[1];
        Assert.assertTrue((first.equals("GTT_AC") && second.equals("GTTGAC")) || (first.equals("GTTGAC") && second.equals("GTT_AC")));
    }
    @Test
    public void swTestHarder1() throws IOException, SubstMatrix.MatrixDimensionException {
        String seq1path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/test2a.txt")).getPath();
        String seq2path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/test2b.txt")).getPath();
        String mxPath = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("matrices/DNAfull.txt")).getPath();
        NWSW swTest = new NWSW();
        List<String[]> result = swTest.smithWaterman(seq1path, seq2path, mxPath, 1);
        String first = result.get(0)[0];
        String second = result.get(0)[1];
        Assert.assertTrue(first.equals("TACGGGCCCGCTA_C") || second.equals("TACGGGCCCGCTA_C"));
        Assert.assertTrue(first.equals("TA___G_CC_CTATC") || second.equals("TA___G_CC_CTATC"));
    }
    @Test
    public void swTestHarder2() throws IOException, SubstMatrix.MatrixDimensionException {
        String seq1path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/test2a.txt")).getPath();
        String seq2path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/test2b.txt")).getPath();
        String mxPath = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("matrices/DNAfull.txt")).getPath();
        NWSW swTest = new NWSW();
        List<String[]> result = swTest.smithWaterman(seq1path, seq2path, mxPath, 1, 5);
        String first = result.get(0)[0];
        String second = result.get(0)[1];
        Assert.assertTrue(first.equals("TACGGGCCCGCTA") || second.equals("TACGGGCCCGCTA"));
        Assert.assertTrue(first.equals("TA___GCC__CTA") || second.equals("TA___GCC__CTA"));
    }
    @Test
    public void nwTest1() throws IOException, SubstMatrix.MatrixDimensionException {
        String seq1path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/testNW1a.txt")).getPath();
        String seq2path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/testNW1b.txt")).getPath();
        String mxPath = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("matrices/nwMXtest.txt")).getPath();
        NWSW nwTest = new NWSW();
        List<String[]> result = nwTest.needlemanWunsch(seq1path, seq2path, mxPath, 4);
        String first = result.get(0)[0];
        String second = result.get(0)[1];
        Assert.assertTrue(first.equals("_ATCGAC") || second.equals("_ATCGAC"));
        Assert.assertTrue(first.equals("CAT__AC") || second.equals("CAT__AC"));
    }
    @Test
    public void nwTest2() throws IOException, SubstMatrix.MatrixDimensionException {
        String seq1path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/testNW2a.txt")).getPath();
        String seq2path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/testNW2b.txt")).getPath();
        String mxPath = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("matrices/nwMCsimple.txt")).getPath();
        NWSW nwTest = new NWSW();
        List<String[]> result = nwTest.needlemanWunsch(seq1path, seq2path, mxPath, 1);
        Assert.assertEquals(3, result.size());
    }
    @Test
    public void nwTest3() throws IOException, SubstMatrix.MatrixDimensionException {
        String seq1path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/testNW2a.txt")).getPath();
        String seq2path = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("sequences/testNW2b.txt")).getPath();
        String mxPath = Objects.requireNonNull(NWSWTest.class.getClassLoader().getResource("matrices/nwMCsimple.txt")).getPath();
        NWSW nwTest = new NWSW();
        List<String[]> result = nwTest.needlemanWunsch(seq1path, seq2path, mxPath, 1);
        Assert.assertTrue(result.get(0)[0].equals("G_ATTACA") || result.get(0)[1].equals("G_ATTACA"));
        Assert.assertTrue((result.get(0)[0].equals(result.get(1)[0]) && result.get(1)[0].equals(result.get(2)[0])) ||
                (result.get(0)[1].equals(result.get(1)[1]) && result.get(1)[1].equals(result.get(2)[1])));
    }
}
