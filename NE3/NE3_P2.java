package paper1.NE0602.NE3;

import org.apache.commons.math3.distribution.TDistribution;

import java.util.ArrayList;
import java.util.Random;

public class NE3_P2 {
    private int k;
    private int d;
    private double alpha;
    private double izpp;
    private ArrayList<Double> izsp = new ArrayList<Double>();  // feasibility tolerance parameter
    private ArrayList<Double> gammas;
    private double feasibleThreshold;
    private ArrayList<Double> ppmu = new ArrayList<Double>();  //mu of the primary performance(standard)
    private ArrayList<Double> ppsigma2 = new ArrayList<Double>();  //sigma2 of primary performance(standard)
    private ArrayList<ArrayList<Double>> spberp = new ArrayList<ArrayList<Double>>();
    private int bestID = -1;
    private double totalSampleSize=0;
    private int goodID = -1;
    private double goodAlternativeSampleSize=0;
    private int goodIDCorrectOrNot=-1; // -1:没找到 0：找到但错了 1：找到且对了

    private int FandG = -1; //1=FG 2=IFG 3=FNG 4=IFNG
    private Random R= new Random();

    public NE3_P2() {
//        R.setSeed(234234);
    }

    public NE3_P2(double alpha, int d, ArrayList<Double> gammas, ArrayList<Double> ppmu, ArrayList<Double> ppsigma2,
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
    public int getFandG() {return FandG;}

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
        int indexgoodalternative = -1;
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
            double[][] sampleSP = new double[I.size()][d];
            for(int ell = 0; ell < n; ell++) {
                for(int i = 0; i < I.size(); i++) {
                    samplePP[i]= R.nextGaussian() * Math.sqrt(ppsigma2.get(I.get(i)))+ppmu.get(I.get(i));
                    for (int j = 0; j < d; j++) {
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

            if(goodID == -1){
                double beta = alpha * (2-Math.PI*Math.PI/6) * (Math.PI*Math.PI/6 - 1);
                double gtildestart = Math.sqrt(-1*t*(Math.log(beta)-Math.log(k-1)-2*Math.log(Math.log(2*t)/Math.log(2)))/2);
                double gstart = gt(t, beta, k-1);
                boolean stepbOrNot = true;
                //update

                outCycle: for(int i = 0;  i < I.size(); i++) {
                    for (int p = 0; p < d; p++) {
                        if (Ztildei[I.get(i)][p] + izsp.get(p)*t < gtildestart) {
                            stepbOrNot = false;
                            break outCycle;

                        }
                    }
                }

                if (stepbOrNot) {
                    boolean findGoodAlternativeOrNot = true;
                    outCycle2: for(int i = 0;  i < I.size(); i++) {
                        for (int j = 0; j < I.size(); j++){
                            if (i!=j){
                                double Sij2 = (Xij2[I.get(i)][I.get(j)] - Xij[I.get(i)][I.get(j)]* Xij[I.get(i)][I.get(j)]/t)/(t-1);
                                if (izpp*t/Math.sqrt(Sij2)-gstart<0){
                                    findGoodAlternativeOrNot = false;
                                    break outCycle2;
                                }
                            }
                        }
                    }
                    if (findGoodAlternativeOrNot){

                        indexgoodalternative = I.get(0);
                        double goodalternativeppp = Xij[I.get(0)][I.get(0)];
                        for (int i = 1;  i < I.size(); i++){
                            if (Xij[I.get(i)][I.get(i)] > goodalternativeppp){
                                goodalternativeppp = Xij[I.get(i)][I.get(i)];
                                indexgoodalternative = I.get(i);
                            }
                        }
                        goodID = indexgoodalternative;
                        goodAlternativeSampleSize = totalSampleSize + I.size()*t;

                        goodIDCorrectOrNot = 1;
                        for (int i = 0;  i < I.size(); i++){
                            boolean fea = true;
                            for (int j = 0; j < d; j++){
                                if (spberp.get(I.get(i)).get(j) < 1 - gammas.get(j)){
                                    fea = false;
                                    break;
                                }
                            }
                            if (ppmu.get(I.get(i)) > ppmu.get(goodID)+izpp && fea){
                                goodIDCorrectOrNot = 0;
                            }
                        }
                        for (int j = 0; j < d; j++){
                            if (spberp.get(goodID).get(j) < 1 - gammas.get(j) - izsp.get(j)){
                                goodIDCorrectOrNot = 0;
                            }
                        }

                    }
                }
            }


            double gtildet = Math.sqrt(-1*t*(Math.log(alpha)-Math.log(k+d-1)-2*Math.log(Math.log(2*t)/Math.log(2)))/2);
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
        if(goodID == -1){
            goodID = bestID;
            goodIDCorrectOrNot = 1;
            goodAlternativeSampleSize = totalSampleSize;
//            System.out.println("已在"+t+"时因找到ba而找到goodalternative，其index为："+goodID+" ss： "+goodAlternativeSampleSize+" goodIDCorrectOrNot: "+goodIDCorrectOrNot);
        }
        long endTime1 = System.currentTimeMillis();
        long runTime1 = endTime1 - startTime1;
        System.out.println("最优方案为："+bestID+" 样本量为："+totalSampleSize+" t: "+t+" 花费时间为："+runTime1);
    }

    public double gt(int t, double beta, int k) {

        double p =  1-beta/((k)*1.0*(Math.log(2*t)/Math.log(2))*(Math.log(2*t)/Math.log(2)));
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
