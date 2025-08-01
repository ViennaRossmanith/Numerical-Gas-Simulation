import java.io.FileWriter;
import java.io.IOException;

public class Simulation {
    private static int nX; // # of grid points in x coordinate direction
    private static int nV; // # of grid points in velocity coordinate direction
    private static double aX; // left endpoint
    private static double bX; // right endpoint
    private static double aV; // min velocity
    private static double bV; // max velocity
    private static double dX; // grid spacing, distance between 2 points on x
    private static double dV; // grid spacing, distance between 2 points on v
    private static double knudsen; // Knudsen #: measures strength of collision operator, closer to zero is
                                   // stronger, closer to infinity is no collisions
    private static double finalTime; // final time
    private static double cfl; // Courant–Friedrichs–Lewy condition
    private static double[] rho; // density in a fluid system
    private static double[] u; // average velocity
    private static double[] T; // temperature
    private static double[] v; // velocity
    private static double[] x; // position
    private static double[][] f; // particle distribution function

    public static void main(String[] args) throws Exception {
        nX = 500;
        nV = 250;
        aX = -1;
        bX = 1;
        aV = -15;
        bV = 15;
        knudsen = 0.1;
        finalTime = 0.32;
        cfl = 1.0;
        dX = (bX - aX) / nX;
        dV = (bV - aV) / nV;
        x = new double[nX];
        v = new double[nV];

        for (int i = 0; i < nX; i++) {
            x[i] = aX + (i + 0.5) * dX;
        }

        for (int i = 0; i < nV; i++) {
            v[i] = aV + (i + 0.5) * dV;
        }

        f = new double[nX][nV];
        rho = new double[nX];
        u = new double[nX];
        T = new double[nX];

        initialCondition();
        outputParams();
        outputMesh();
        outputF(0);
        computeMoments();
        outputMoments(0);
        timeAdvance();
        outputMoments(1);
        outputF(1);
    }

    public static void initialCondition() {
        double p_rho;
        double p_u;
        double p_T;
        for (int i = 0; i < nX; i++) {
            for (int j = 0; j < nV; j++) {
                if (x[i] < 0) {
                    p_rho = 2.25;
                    p_u = 0.0;
                    p_T = 1.125;
                    f[i][j] = p_rho / (Math.sqrt(2 * Math.PI * p_T)) * Math.exp(-Math.pow(v[j] - p_u, 2) / (2 * p_T));
                } else {
                    p_rho = 3.0 / 7.0;
                    p_u = 0.0;
                    p_T = 1.0 / 6.0;
                    f[i][j] = p_rho / (Math.sqrt(2 * Math.PI * p_T)) * Math.exp(-Math.pow(v[j] - p_u, 2) / (2 * p_T));
                }
            }
        }
    }

    public static void outputParams() throws IOException {
        FileWriter paramWriter = new FileWriter("Parameters.txt");
        paramWriter.write("nX: " + nX + "\n");
        paramWriter.write("nV: " + nV + "\n");
        paramWriter.write("aX: " + String.format("%.15e", aX) + "\n");
        paramWriter.write("bX: " + String.format("%.15e", bX) + "\n");
        paramWriter.write("aV: " + String.format("%.15e", aV) + "\n");
        paramWriter.write("bV: " + String.format("%.15e", bV) + "\n");
        paramWriter.write("Knudsen Number: " + String.format("%.15e", knudsen) + "\n");
        paramWriter.write("Final Time: " + String.format("%.15e", finalTime));
        paramWriter.close();
    }

    public static void outputMesh() throws IOException {
        FileWriter xWriter = new FileWriter("OutputX.txt");
        FileWriter vWriter = new FileWriter("OutputV.txt");

        for (int i = 0; i < nX; i++) {
            xWriter.write(String.format("%.15e", x[i]) + "\n");
        }
        xWriter.close();

        for (int i = 0; i < nV; i++) {
            vWriter.write(String.format("%.15e", v[i]) + "\n");
        }
        vWriter.close();
    }

    public static void outputF(int timeIndex) throws IOException {
        String fileName = "OutputF" + String.valueOf(timeIndex) + ".txt";
        FileWriter fWriter = new FileWriter(fileName);
        for (int i = 0; i < nX; i++) {
            for (int j = 0; j < nV; j++) {
                fWriter.write(String.format("%.15e\n", f[i][j]));
            }
        }
        fWriter.close();
    }

    public static void outputMoments(int timeIndex) throws IOException {
        String fileNameRho = "OutputRho" + String.valueOf(timeIndex) + ".txt";
        String fileNameU = "OutputU" + String.valueOf(timeIndex) + ".txt";
        String fileNameT = "OutputT" + String.valueOf(timeIndex) + ".txt";

        FileWriter outputRho = new FileWriter(fileNameRho);
        FileWriter outputU = new FileWriter(fileNameU);
        FileWriter outputT = new FileWriter(fileNameT);
        for (int i = 0; i < rho.length; i++) {
            outputRho.write(String.format("%.15e\n", rho[i]));
        }
        outputRho.close();

        for (int i = 0; i < u.length; i++) {
            outputU.write(String.format("%.15e\n", u[i]));
        }
        outputU.close();

        for (int i = 0; i < T.length; i++) {
            outputT.write(String.format("%.15e\n", T[i]));
        }
        outputT.close();
    }

    public static void computeMoments() {
        for (int i = 0; i < nX; i++) {
            rho[i] = 0;
            for (int j = 0; j < nV; j++) {
                rho[i] += dV * f[i][j];
            }
        }
        for (int i = 0; i < nX; i++) {
            u[i] = 0;
            for (int j = 0; j < nV; j++) {
                u[i] += dV * v[j] * f[i][j];
            }
            u[i] = u[i] / rho[i];
        }
        for (int i = 0; i < nX; i++) {
            T[i] = 0;
            for (int j = 0; j < nV; j++) {
                T[i] += dV * Math.pow(v[j] - u[i], 2) * f[i][j];
            }
            T[i] = T[i] / rho[i];
        }
    }

    public static void timeAdvance() {
        double vMax = Math.max(Math.abs(aV), Math.abs(bV));
        // double vMin = v[nV/2];
        double dT = dX * cfl / vMax;
        int nSteps = (int) Math.ceil(finalTime / dT);
        dT = finalTime / nSteps;
        int[] m = new int[nV];
        double[] nuBar = new double[nV];

        for (int j = 0; j < nV; j++) {
            // System.out.println("j:" + j + " nu:" + dT*v[j]/dX);
            double nu = dT * Math.abs(v[j]) / dX;
            m[j] = (int) Math.floor(nu);
            nuBar[j] = nu - m[j];
            System.out.println(nuBar[j]);
        }

        double[][] fOld = f;
        for (int n = 1; n <= nSteps; n++) {
            for (int i = 0; i < nX; i++) {
                for (int j = 0; j < nV / 2; j++) {
                    int shift = i + m[j];
                    int shiftP1 = shift + 1;
                    if (shift > nX - 1) {
                        shift = nX - 1;
                    }
                    if (shiftP1 > nX - 1) {
                        shiftP1 = nX - 1;
                    }
                    f[i][j] = fOld[shift][j] - nuBar[j] * (fOld[shift][j] - fOld[shiftP1][j]);
                }
            }
            for (int i = 0; i < nX; i++) {
                for (int j = nV / 2; j < nV; j++) {
                    int shift = i - m[j];
                    int shiftM1 = shift - 1;
                    if (shift < 0) {
                        shift = 0;
                    }
                    if (shiftM1 < 0) {
                        shiftM1 = 0;
                    }
                    f[i][j] = fOld[shift][j] - nuBar[j] * (fOld[shift][j] - fOld[shiftM1][j]);
                }
            }
            /*
             * double[][] fOld = f;
             * for (int n = 1; n <= nSteps; n++) {
             * for (int i = 0; i < nX; i++) {
             * for (int j = 0; j < nV; j++) {
             * int shift = i + m[j];
             * int shiftP1 = shift + 1;
             * int shiftM1 = shift - 1;
             * if (shift > nX - 1) {
             * shift = nX - 1;
             * }
             * if (shiftP1 > nX - 1) {
             * shiftP1 = nX - 1;
             * }
             * if (shiftM1 < 0) {
             * shiftM1 = 0;
             * }
             * double nu = dT * v[j] / dX;
             * f[i][j] = fOld[shift][j] - 0.5 * nu * (fOld[shiftP1][j] - fOld[shiftM1][j])
             * + 0.5 * nu * nu * (fOld[shiftP1][j] - 2.0 * fOld[shift][j] +
             * fOld[shiftM1][j]);
             * }
             * }
             */
            computeMoments();
            for (int i = 0; i < nX; i++) {
                for (int j = 0; j < nV; j++) {
                    double MM = rho[i] / Math.sqrt(2 * Math.PI * T[i])
                            * Math.exp(-(Math.pow(v[j] - u[i], 2) / (2 * T[i])));
                    f[i][j] = (knudsen / (knudsen + dT)) * f[i][j] + (dT / (knudsen + dT)) * MM;
                }
            }
            fOld = f;
        }
    }
}
