package paper1.NE0602.NE2;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;

import java.util.ArrayList;
import java.util.Random;

public class NE2_P2 {
    private int k;
    private int d;
    private double alpha;
    private double izpp;
    private ArrayList<Double> izsp = new ArrayList<Double>();  // feasibility tolerance parameter
    private ArrayList<Double> gammas;
    private ArrayList<Double> ppmu = new ArrayList<Double>();  //mu of the primary performance(normal)
    private ArrayList<Double> ppsigma2 = new ArrayList<Double>();  //sigma2 of primary performance(normal)
    private ArrayList<ArrayList<Double>> spberp = new ArrayList<ArrayList<Double>>();
    private int bestID = -1;
    private double totalSampleSize=0;
    private int goodID = -1;
    private double goodAlternativeSampleSize=0;
    private int goodIDCorrectOrNot=-1; // -1:没找到 0：找到但错了 1：找到且对了

    private Random R= new Random();

    public NE2_P2() {
//        R.setSeed(234234);
    }

    public NE2_P2(double alpha, int d, ArrayList<Double> gammas, ArrayList<Double> ppmu, ArrayList<Double> ppsigma2,
                  ArrayList<ArrayList<Double>> spberp, long seed, double izpp, ArrayList<Double> izsp){
        this.alpha = alpha;
        this.d = d;
        this.gammas = gammas;
        this.ppmu.addAll(ppmu);
        this.ppsigma2.addAll(ppsigma2);
        this.spberp.addAll(spberp);
        R.setSeed(seed);
        this.izpp = izpp;
        this.izsp.addAll(izsp);
    }

    public int getBestID() {
        return bestID;
    }
    public int getGoodID() {
        return goodID;
    }

    public double getTotalSampleSize() {
        return totalSampleSize;
    }
    public double getGoodAlternativeSampleSize() {
        return goodAlternativeSampleSize;
    }
    public int getGoodIDCorrectOrNot() {
        return goodIDCorrectOrNot;
    }

    public void run() {
        long startTime1 = System.currentTimeMillis();
        k = ppmu.size();
        ArrayList<Integer> I = new ArrayList<Integer>();
        for(int i = 0 ; i < k ; i++) {
            I.add(i);
        }

        double[][] Xij =  new double[k][k];
        double[][] Xij2 = new double[k][k];
        double[][] Ztildei = new double[k][d];

        int t = 1;
        int r = 1;
        int n = 0;
        while(I.size()>1) {
            System.out.println("t "+t+" I.SIZE "+I.size());
            if(t == 1) {
                t = t + 1;
                n = 2;
            }else{
                n = t;
                r = r + 1;
                t = 2*t;
            }
            double[] samplePP = new double[I.size()];
            double[][] sampleSP = new double[I.size()][d];
            for(int ell = 0; ell < n; ell++) {
                for(int i = 0; i < I.size(); i++) {
                    //primary performance measure sample
                    samplePP[i] = R.nextGaussian() * Math.sqrt(ppsigma2.get(I.get(i)))+ppmu.get(I.get(i));
                    //secondary performance measure sample
                    for(int j = 0; j < d; j++) {
                        sampleSP[i][j] = bernoullisample(spberp.get(I.get(i)).get(j));
                    }
                }

                for(int i = 0; i < I.size(); i++) {
                    Xij[I.get(i)][I.get(i)] = Xij[I.get(i)][I.get(i)] + samplePP[i];
                    Xij2[I.get(i)][I.get(i)] = Xij2[I.get(i)][I.get(i)] + samplePP[i] * samplePP[i];

                    for(int j = i + 1; j < I.size(); j++) {
                        Xij[I.get(i)][I.get(j)] = Xij[I.get(i)][I.get(j)]+samplePP[i]-samplePP[j];
                        Xij[I.get(j)][I.get(i)] = Xij[I.get(j)][I.get(i)]+samplePP[j]-samplePP[i];
                        Xij2[I.get(i)][I.get(j)] = Xij2[I.get(i)][I.get(j)] +(samplePP[i]-samplePP[j])*(samplePP[i]-samplePP[j]);
                        Xij2[I.get(j)][I.get(i)] = Xij2[I.get(i)][I.get(j)];
                    }

                    for (int k = 0; k < d; k++){
                        Ztildei[I.get(i)][k] = Ztildei[I.get(i)][k] + sampleSP[i][k] - 1 + gammas.get(k);
                    }
                }
            }


            double gtildet = Math.sqrt(-1.0*t*(Math.log(alpha)-Math.log(k+d-1)-2.0*Math.log(Math.log(2.0*t)/Math.log(2.0)))/2);
            double gt = gt(t, alpha, k+d-1);
            int[] feasibleOrNot = new int[k]; // 1=infeasible 0=feasible
            for(int i = 0;  i < I.size(); i++) {
                for (int p = 0; p < d; p++) {
                    if (Ztildei[I.get(i)][p] <= -gtildet) {
                        feasibleOrNot[I.get(i)] = 1;
                        break;
                    } else if (Ztildei[I.get(i)][p] < gtildet) {
                        feasibleOrNot[I.get(i)] = -1;
                    }
                }
            }
            for(int i = 0;  i < I.size(); i++) {
                if (feasibleOrNot[I.get(i)] == 1){
                    I.remove(i);
                    i--;
                    totalSampleSize = totalSampleSize + t;
                }else{
                    for (int j = 0; j < I.size(); j++) {
                        double Sij2 = (Xij2[I.get(i)][I.get(j)] - Xij[I.get(i)][I.get(j)]* Xij[I.get(i)][I.get(j)]/t)/(t-1);
                        double Zji = Xij[I.get(j)][I.get(i)]/Math.sqrt(Sij2);
                        if(i!=j) {
                            if (feasibleOrNot[I.get(j)] == 0 && Zji > gt){
                                I.remove(i);
                                i--;
                                totalSampleSize = totalSampleSize + t;
                                break;
                            }
                        }
                    }
                }
            }
            if(I.size()==1) {
                totalSampleSize = totalSampleSize + t;
            }

        }
        bestID = I.get(0);
        long endTime1 = System.currentTimeMillis();
        long runTime1 = endTime1 - startTime1;
    }

    public double gt(long t, double beta, int k) {

//        System.out.println(">>> DEBUG: t = " + t + ", beta = " + beta + ", k = " + k);
        if (t <= 0 || k <= 0 || Double.isNaN(beta)) {
            throw new IllegalArgumentException("t, k must be positive and beta must be a number.");
        }

        double log2_2t = Math.log(2.0 * t) / Math.log(2);
//        System.out.println(">>> log2_2t = " + log2_2t);

        if (Double.isInfinite(log2_2t) || Double.isNaN(log2_2t)) {
            throw new IllegalStateException("log2_2t is not finite");
        }

        double p = 1 - beta / (k * log2_2t * log2_2t);
//        System.out.println(">>> p = " + p);

        if (Double.isNaN(p) || p < 0 || p > 1) {
            throw new IllegalArgumentException("Invalid p: " + p);
        }
//        System.out.println(p);

//        double p2 =  1-beta/((k)*1.0*(Math.log(2.0*t)/Math.log(2.0))*(Math.log(2.0*t)/Math.log(2)));
//        System.out.println(p+" "+p2);
        TDistribution td = new TDistribution(t-1);
//        System.out.println(t+" "+p2+" "+beta);
        double  gt = Math.sqrt(t)*td.inverseCumulativeProbability(p);

        return gt;
    }

    public double bernoullisample(double p){
//        System.out.println(R.nextDouble());
        boolean result = R.nextDouble() < p;

        if (result) {
            return 1;
        } else {
            return 0;
        }
    }

    public double exponentialsample(double lambda){
        double u = R.nextDouble();
        double exponentialRandom = -Math.log(1 - u) / lambda;
        return exponentialRandom;
    }

    public static double inverseCumulativeProbabilityforSTD(double p){
        NormalDistribution normalDistribution = new NormalDistribution();
        return normalDistribution.inverseCumulativeProbability(p);
    }
}
