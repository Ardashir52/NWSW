package nprg013.zapoctovy_program.nwsw;

import java.util.List;

/**
 * Auxiliary class for gap value computations.
 */
public class GapPenalty {
    /**
     * Contain one or two doubles according to the type of gap penalty.
     */
    private final double[] gapPenalty;

    /**
     * Initializes the class according to the number of numbers passed to it.
     *
     * @param gapValue  The list of values used to calculate gap penalties.
     * @throws GapPenaltyNumericalError  This is thrown in case of input parsing failure.
     */
    public GapPenalty(List<String> gapValue) {
        try {
            if (gapValue.size() > 1) {
                this.gapPenalty = new double[]{Double.parseDouble(gapValue.get(0)), Double.parseDouble(gapValue.get(1))};
            } else if (!gapValue.isEmpty()) {
                this.gapPenalty = new double[]{Double.parseDouble(gapValue.get(0))};
            } else {
                this.gapPenalty = new double[]{3.0};
            }
        }
        catch (NumberFormatException ne) {
            throw new GapPenaltyNumericalError("Wrong gap penalty number format!");
        }
    }

    /**
     * The constructor for use in interactive mode.
     *
     * @param gapValue  List is not needed as only one value is accepted for simplicity.
     * @throws GapPenaltyNumericalError  This is thrown in case of parsing failure.
     */
    public GapPenalty(String gapValue) {
        try {
            this.gapPenalty = new double[] {Double.parseDouble(gapValue)};
        }
        catch (NumberFormatException ne) {
            throw new GapPenaltyNumericalError("Wrong gap penalty number format!");
        }
    }
    /**
     *
     * @return  The value of the linear part of gap penalty.
     */
    public double linearPart() {
        return gapPenalty[0];
    }

    /**
     * Simple method to count gap value according to the node distance.
     *
     * @param k  The distance between nodes parameter
     * @return  The value of gap stretching k nodes away.
     */
    public double countAffine(int k) {
        return (gapPenalty[1] + k * gapPenalty[0]);
    }

    /**
     * The method used to determine whether linear or affine gap should be used in calculations.
     *
     * @return  The type of gap.
     */
    public Type getType() {
        switch (this.gapPenalty.length) {
            case 1:
                return Type.LINEAR;
            case 2:
                return Type.AFFINE;
            default:
                return Type.ERROR;
        }
    }

    /**
     * Enumerates the possible types of a gap penalty.
     */
    public enum Type {
        LINEAR,
        AFFINE,
        ERROR
    }

    /**
     * Exception raised in the case that values provided for gap value cannot be parsed to double.
     */
    public static class GapPenaltyNumericalError extends NumberFormatException {
        public GapPenaltyNumericalError(String s) {
            super(s);
        }
    }
}
