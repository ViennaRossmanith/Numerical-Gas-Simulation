import java.io.FileWriter;
import java.io.IOException;

public class AlternateCollision {
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

    public static void main(String[] args) throws Exception {
        nX = 256;
        nV = 128;
        knudsen = 0.01;
        finalTime = 0.16;
        cfl = 1.95;
        aX = -1.25;
        bX = 1.25;
        aV = -7.0;
        bV = 7.0;
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

        im2v[0] = nX - 2;
        im2v[1] = nX - 1;
        im1v[0] = nX - 1;
        ip1v[nX - 1] = 0;
        ip2v[nX - 2] = 0;
        ip2v[nX - 1] = 1;

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
                if (Math.abs(x[i]) < 0.5) {
                    p_rho = 1.000;
                    p_u = 0.250;
                    p_T = 1.000;
                    f[i][j] = p_rho / (Math.sqrt(2 * Math.PI * p_T)) * Math.exp(-Math.pow(v[j] - p_u, 2) / (2 * p_T));
                } else {
                    p_rho = 0.125;
                    p_u = -0.10;
                    p_T = 0.800;
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
        for (int i = 0; i < nX; i++) {
            int im2 = im2v[i];
            int im1 = im1v[i];
            int ip1 = ip1v[i];
            int ip2 = ip2v[i];
            for (int j = 0; j < nV; j++) {
                double nuBar = (v[j] * p_dT) / dX;
                if (nuBar > 0) {
                    f[i][j] = fOld[i][j] - (nuBar / 6)
                            * (fOld[im2][j] - 6 * fOld[im1][j] + 3 * fOld[i][j] + 2 * fOld[ip1][j])
                            + (Math.pow(nuBar, 2) / 2) * (fOld[im1][j] - 2 * fOld[i][j] + fOld[ip1][j])
                            - (Math.pow(nuBar, 3) / 6)
                                    * (-fOld[im2][j] + 3 * fOld[im1][j] - 3 * fOld[i][j] + fOld[ip1][j]);
                } else if (nuBar < 0) {
                    f[i][j] = fOld[i][j] - (nuBar / 6)
                            * (-2 * fOld[im1][j] - 3 * fOld[i][j] + 6 * fOld[ip1][j] - fOld[ip2][j])
                            + (Math.pow(nuBar, 2) / 2) * (fOld[im1][j] - 2 * fOld[i][j] + fOld[ip1][j])
                            - (Math.pow(nuBar, 3) / 6)
                                    * (-fOld[im1][j] + 3 * fOld[i][j] - 3 * fOld[ip1][j] + fOld[ip2][j]);
                }
            }
        }
    }

    public static void collision(double dT) {
        double[][] MM = new double[nX][nV];
        double mu;
        double theta = (dT * (dT + 12 * knudsen)) / ((dT + 3 * knudsen) * (dT + 4 * knudsen));

        for (int i = 0; i < nX; i++) {
            double factor = (dV / Math.sqrt(2 * Math.PI * T[i]));
            double A0 = 0;
            double A1 = 0;
            double A2 = 0;
            double A3 = 0;
            double A4 = 0;
            for (int j = 0; j < nV; j++) {
                mu = (v[j] - u[i]) / Math.sqrt(T[i]);
                A0 += Math.pow(mu, 0) * Math.exp(-0.5 * (Math.pow(mu, 2)));
                A1 += Math.pow(mu, 1) * Math.exp(-0.5 * (Math.pow(mu, 2)));
                A2 += Math.pow(mu, 2) * Math.exp(-0.5 * (Math.pow(mu, 2)));
                A3 += Math.pow(mu, 3) * Math.exp(-0.5 * (Math.pow(mu, 2)));
                A4 += Math.pow(mu, 4) * Math.exp(-0.5 * (Math.pow(mu, 2)));
            }
            A0 *= factor;
            A1 *= factor;
            A2 *= factor;
            A3 *= factor;
            A4 *= factor;
            double d = Math.pow(A2, 3) - 2 * A1 * A2 * A3 + A0 * Math.pow(A3, 2) + Math.pow(A1, 2) * A4 - A0 * A2 * A4;
            double[] a = new double[3];
            a[0] = (1 / d) * (Math.pow(A1, 2) + A2 * (2 * A2 - A0 - A4) - A3 * (2 * A1 - A3));
            a[1] = (1 / d) * (A1 * (A4 - A2) + A3 * (A0 - A2));
            a[2] = (1 / d) * (A1 * (A1 - A3) + A2 * (A2 - A0));

            for (int j = 0; j < nV; j++) {
                mu = (v[j] - u[i]) / Math.sqrt(T[i]);
                MM[i][j] = (rho[i] / Math.sqrt(2 * Math.PI * T[i])) * Math.exp(-(Math.pow(mu, 2) / 2))
                        * (a[0] + a[1] * mu + a[2] * (Math.pow(mu, 2) - 1));
                f[i][j] = theta * MM[i][j] + (1 - theta) * f[i][j];
            }
        }
    }

    public static void conservation(double[] totals) {
        totals[0] = 0;
        totals[1] = 0;
        totals[2] = 0;
        for (int i = 0; i < nX; i++) {
            for (int j = 0; j < nV; j++) {
                totals[0] += dX * dV * f[i][j];
                totals[1] += dX * dV * v[j] * f[i][j];
                totals[2] += dX * dV * Math.pow(v[j], 2) * f[i][j];
            }
        }
    }

    public static void outputConservation(double[] time, double[] mass, double[] momentum, double[] energy, int length)
            throws IOException {
        String fileName = "OutputConservation.txt";

        FileWriter outputConservation = new FileWriter(fileName);

        for (int i = 0; i < length; i++) {
            outputConservation.write(
                    String.format("%.15e    %.15e    %.15e    %.15e%n",
                            time[i], mass[i], momentum[i], energy[i]));
        }

        outputConservation.close();
    }

    public static void timeAdvance() throws IOException {
        double[][] fOld = new double[nX][nV];
        double vMax = Math.max(Math.abs(aV), Math.abs(bV));
        double dT = dX * cfl / vMax;
        int nSteps = (int) Math.ceil(finalTime / dT);
        dT = finalTime / nSteps;
        cfl = vMax * dT / dX;
        double[] time = new double[nSteps + 1];
        double[] mass = new double[nSteps + 1];
        double[] mom = new double[nSteps + 1];
        double[] nrg = new double[nSteps + 1];
        double[] totals = new double[3];
        time[0] = 0;
        conservation(totals);
        mass[0] = totals[0];
        mom[0] = totals[1];
        nrg[0] = totals[2];

        for (int n = 1; n <= nSteps; n++) {
            computeMoments();
            collision(dT / 2);
            for (int i = 0; i < nX; i++) {
                for (int j = 0; j < nV; j++) {
                    fOld[i][j] = f[i][j];
                }
            }
            transport(dT, fOld);
            computeMoments();
            collision(dT / 2);
            time[n] = n * dT;
            conservation(totals);
            mass[n] = totals[0];
            mom[n] = totals[1];
            nrg[n] = totals[2];
        }

        String fileName = "Steps.txt";
        FileWriter nStepsWriter = new FileWriter(fileName);
        nStepsWriter.write(String.valueOf(nSteps));
        nStepsWriter.close();

        outputConservation(time, mass, mom, nrg, nSteps + 1);
    }
}
