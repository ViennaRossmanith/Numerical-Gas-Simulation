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
    private static int[] im2v; // spatial indeces
    private static int[] im1v;
    private static int[] ip1v;
    private static int[] ip2v;
    private static double[][] fOld;

    public static void main(String[] args) throws Exception {
        nX = 500;
        nV = 250;
        aX = -1;
        bX = 1;
        aV = -15;
        bV = 15;
        knudsen = 0.001;
        finalTime = 0.32;
        cfl = 1.0;
        dX = (bX - aX) / nX;
        dV = (bV - aV) / nV;
        x = new double[nX];
        v = new double[nV];
        im2v = new int[nX];
        im1v = new int[nX];
        ip1v = new int[nX];
        ip2v = new int[nX];

        for (int i = 0; i < nX; i++) {
            im2v[i] = i - 2;
            im1v[i] = i - 1;
            ip1v[i] = i + 1;
            ip2v[i] = i + 2;
        }

        im2v[0] = 0;
        im2v[1] = 0;
        im1v[0] = 0;
        ip1v[nX - 1] = nX - 1;
        ip2v[nX - 2] = nX - 1;
        ip2v[nX - 1] = nX - 1;

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

    public static void transport(double p_dT, double[][] fOld) {
        double[] nuBar = new double[nV];
        for (int i = 0; i < nX; i++) {
            for (int j = 0; j < nV; j++) {
                nuBar[j] = (v[j] * p_dT) / dX;
                if (nuBar[j] > 0) {
                    f[i][j] = fOld[i][j] - (nuBar[j] / 6)
                            * (fOld[im2v[i]][j] - 6 * fOld[im1v[i]][j] + 3 * (fOld[i][j]) + 2 * fOld[ip1v[i]][j])
                            + (Math.pow(nuBar[j], 2) / 2) * (fOld[im1v[i]][j] - 2 * fOld[i][j] + fOld[ip1v[i]][j])
                            - (Math.pow(nuBar[j], 3) / 6)
                                    * (-fOld[im2v[i]][j] + 3 * fOld[im1v[i]][j] - 3 * fOld[i][j] + fOld[ip1v[i]][j]);
                } else if (nuBar[j] < 0) {
                    f[i][j] = fOld[i][j] - (nuBar[j] / 6)
                            * (-2 * fOld[im1v[i]][j] - 3 * fOld[i][j] + 6 * (fOld[ip1v[i]][j]) - fOld[ip2v[i]][j])
                            + (Math.pow(nuBar[j], 2) / 2) * (fOld[im1v[i]][j] - 2 * fOld[i][j] + fOld[ip1v[i]][j])
                            - (Math.pow(nuBar[j], 3) / 6)
                                    * (-fOld[im1v[i]][j] + 3 * fOld[i][j] - 3 * fOld[ip1v[i]][j] + fOld[ip2v[i]][j]);
                }
            }
        }
    }

    public static void collision(double dT) {
        double[][] MM = new double[nX][nV];
        double[][] mu = new double[nX][nV];
        double theta = (dT * (dT + 12 * knudsen)) / ((dT + 3 * knudsen) * (dT + 4 * knudsen));
        for (int i = 0; i < nX; i++) {
            for (int j = 0; j < nV; j++) {
                mu[i][j] = (v[j] - u[i]) / Math.sqrt(T[i]);
                MM[i][j] = (rho[i] / Math.sqrt(2 * Math.PI * T[i])) * Math.exp(-(Math.pow(mu[i][j], 2) / 2));
                f[i][j] = theta * MM[i][j] + (1 - theta) * fOld[i][j];
            }
        }
    }

    public static void timeAdvance() {
        double vMax = Math.max(Math.abs(aV), Math.abs(bV));
        double dT = dX * cfl / vMax;
        int nSteps = (int) Math.ceil(finalTime / dT);
        dT = finalTime / nSteps;

        fOld = f;
        transport(dT / 2, fOld);
        computeMoments();
        collision(dT);
        fOld = f;
        for (int n = 2; n <= nSteps; n++) {
            transport(dT, fOld);
            computeMoments();
            collision(dT);
            fOld = f;
        }
        transport(dT / 2, fOld);
        computeMoments();
    }
}
