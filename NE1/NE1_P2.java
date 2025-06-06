package paper1.NE0602.NE1;

import org.apache.commons.math3.distribution.TDistribution;

import java.util.ArrayList;
import java.util.Random;

public class NE1_P2 {
    private int k;
    private double alpha;
    private double gamma;
    private double feasibleThreshold;
    private ArrayList<Double> ppmu = new ArrayList<Double>();  //mu of the primary performance(standard)
    private ArrayList<Double> ppsigma2 = new ArrayList<Double>();  //sigma2 of primary performance(standard)
    private ArrayList<Double> spberp = new ArrayList<Double>();
    private int bestID = -1;
    private double totalSampleSize=0;

    private Random R= new Random();

    public NE1_P2() {
//        R.setSeed(234234);
    }

    public NE1_P2(double alpha, double gamma, double feasibleThreshold, ArrayList<Double> ppmu, ArrayList<Double> ppsigma2,
                  ArrayList<Double> spberp, long seed){
        this.alpha = alpha;
        this.gamma = gamma;
        this.feasibleThreshold = feasibleThreshold;
        this.ppmu.addAll(ppmu);
        this.ppsigma2.addAll(ppsigma2);
        this.spberp.addAll(spberp);
        R.setSeed(seed);
    }

    public int getBestID() {
        return bestID;
    }

    public double getTotalSampleSize() {
        return totalSampleSize;
    }

    public void run() {
        long startTime1 = System.currentTimeMillis();
        k  = ppmu.size();
        ArrayList<Integer> I = new ArrayList<Integer>();
        for(int i = 0 ; i < k ; i++) {
            I.add(i);
        }

        double[][] Xij =  new double[k][k];
        double[][] Xij2 = new double[k][k];
        double[] Ztildei = new double[k];

        int t = 1;
        int r = 1;
        int n = 0;
        while(I.size()>1) {
            if(t == 1) {
                t = t + 1;
                n = 2;
            }else{
                n = t;
                r = r + 1;
                t = 2*t;
            }
            double[] samplePP = new double[I.size()];
            double[] sampleSP = new double[I.size()];
            for(int ell = 0; ell < n; ell++) {
                for(int i = 0; i < I.size(); i++) {
                    samplePP[i]= R.nextGaussian() * Math.sqrt(ppsigma2.get(I.get(i)))+ppmu.get(I.get(i));
                    sampleSP[i]= bernoullisample(spberp.get(I.get(i)));
                }

                for(int i = 0; i < I.size(); i++) {
                    Xij[I.get(i)][I.get(i)] = Xij[I.get(i)][I.get(i)] + samplePP[i];
                    Xij2[I.get(i)][I.get(i)] = Xij2[I.get(i)][I.get(i)] + samplePP[i] * samplePP[i];
                    Ztildei[I.get(i)] = Ztildei[I.get(i)] + sampleSP[i] - 1 + gamma;

                    for(int j = i + 1; j < I.size(); j++) {
                        Xij[I.get(i)][I.get(j)] = Xij[I.get(i)][I.get(j)]+samplePP[i]-samplePP[j];
                        Xij[I.get(j)][I.get(i)] = Xij[I.get(j)][I.get(i)]+samplePP[j]-samplePP[i];
                        Xij2[I.get(i)][I.get(j)] = Xij2[I.get(i)][I.get(j)] +(samplePP[i]-samplePP[j])*(samplePP[i]-samplePP[j]);
                        Xij2[I.get(j)][I.get(i)] = Xij2[I.get(i)][I.get(j)];
                    }
                }
            }

            double gtildet = Math.sqrt(-1.0*t*(Math.log(alpha)-Math.log(k)-2.0*Math.log(Math.log(2.0*t)/Math.log(2.0)))/2.0);
            double gt = gt(t);
            for(int i = 0;  i < I.size(); i++) {
                if (Ztildei[I.get(i)] <= -gtildet){
                    I.remove(i);
                    i--;
                    totalSampleSize = totalSampleSize + t;
                }else{
                    for (int j = 0; j < I.size(); j++) {
                        double Sij2 = (Xij2[I.get(i)][I.get(j)] - Xij[I.get(i)][I.get(j)]* Xij[I.get(i)][I.get(j)]/t)/(t-1);
                        double Zji = Xij[I.get(j)][I.get(i)]/Math.sqrt(Sij2);
                        if(i!=j) {
                            if (Ztildei[I.get(j)] >= gtildet && Zji > gt){
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

    public double gt(int t) {

        double p =  1-alpha/(k*1.0*(Math.log(2.0*t)/Math.log(2.0))*(Math.log(2.0*t)/Math.log(2)));
        TDistribution td = new TDistribution(t-1);

        double  gt = Math.sqrt(t)*td.inverseCumulativeProbability(p);

        return gt;
    }

    public double bernoullisample(double p){
        boolean result = R.nextDouble() < p;

        if (result) {
            return 1;
        } else {
            return 0;
        }
    }
}
