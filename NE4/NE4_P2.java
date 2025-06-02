package paper1.NE0602.NE4;

import org.apache.commons.math3.distribution.TDistribution;

import java.util.ArrayList;
import java.util.Random;

public class NE4_P2 {
    private int k;
    private int d;
    private double alpha;
    private ArrayList<Double> gammas;
    private double feasibleThreshold;
    private ArrayList<int[]> listALT;
    private int bestID = -1;
    private double totalSampleSize=0;

    private Random R= new Random();

    public NE4_P2() {
//        R.setSeed(234234);
    }

    public NE4_P2(double alpha, int d, ArrayList<Double> gammas, ArrayList<int[]> listALT, long seed){
        this.alpha = alpha;
        this.d = d;
        this.gammas = gammas;
        this.listALT = listALT;
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
        k = listALT.size();
        NE4_ED simulation = new NE4_ED(R.nextLong());
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
        int sampleSizeFC = 0;
        int sampleSizeOC = 0;
        int fcEliminate = 0;
        while(I.size()>1) {
//            System.out.println(t);
            if(t == 1) {
                t = t + 1;
                n = 2;
            }else{
                n = t;
                r = r + 1;
                t = 2*t;
            }
            System.out.println(t+" I.size "+I.size());
            double[] samplePP = new double[I.size()];
            double[][] sampleSP = new double[I.size()][d];
            for(int ell = 0; ell < n; ell++) {
                for(int i = 0; i < I.size(); i++) {

                    simulation.setAlt(listALT.get(I.get(i)));
                    simulation.run();
                    samplePP[i]= -simulation.getAverageWTForC4C5();
                    sampleSP[i][0] = simulation.getSatisfyC1OrNot();
                    sampleSP[i][1] = simulation.getSatisfyC2OrNot();
                    sampleSP[i][2] = simulation.getSatisfyC3OrNot();
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

            double gtildet = Math.sqrt(-1*t*(Math.log(alpha)-Math.log(k+d-1)-2*Math.log(Math.log(2*t)/Math.log(2)))/2);
            double gt = gt(t);
            int[] feasibleOrNot = new int[k]; // 1=infeasible 0=feasible
            for(int i = 0;  i < I.size(); i++) {
                for (int p = 0; p < d; p++) {
                    if (Ztildei[I.get(i)][p] <= -gtildet) {
                        feasibleOrNot[I.get(i)] = 1;
//                        System.out.println(t+"时方案"+I.get(i)+"因约束"+p+"表现为"+Ztildei[I.get(i)][p]/t+"小于"+-gtildet/t+"被判定为不可行");
                        break;
                    } else if (Ztildei[I.get(i)][p] < gtildet) {
                        feasibleOrNot[I.get(i)] = -1;
                    }
                }
            }
            for(int i = 0;  i < I.size(); i++) {
                if (feasibleOrNot[I.get(i)] == 1){
//                    fcEliminate=fcEliminate+1;
//                    System.out.println(t+"时通过可行性分析去掉了 "+I.get(i)+" 可行性分析已去掉 "+fcEliminate);
                    I.remove(i);
                    i--;
                    totalSampleSize = totalSampleSize + t;
                    sampleSizeFC +=t;
                }else{
                    for (int j = 0; j < I.size(); j++) {
//                        System.out.println(I.get(i)+" "+I.get(j));
//                        System.out.println(i+" "+j);
                        double Sij2 = (Xij2[I.get(i)][I.get(j)] - Xij[I.get(i)][I.get(j)]* Xij[I.get(i)][I.get(j)]/t)/(t-1);
                        double Zji = Xij[I.get(j)][I.get(i)]/Math.sqrt(Sij2);
                        if(i!=j) {
                            if (feasibleOrNot[I.get(j)] == 0 && Zji > gt){
//                                System.out.println(t+"时"+I.get(j)+"凭借pp为"+Xij[I.get(j)][I.get(j)]/t+"的表现"+"通过最优性分析去掉了"+I.get(i)
//                                        +" 其pp为："+Xij[I.get(i)][I.get(i)]/t);
                                I.remove(i);
                                i--;
                                totalSampleSize = totalSampleSize + t;
                                sampleSizeOC +=t;
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
        System.out.println("最优方案为："+bestID+" 样本量为："+totalSampleSize+" 花费时间为："+runTime1);
        System.out.println("可行性分析花费"+sampleSizeFC+" 占比"+sampleSizeFC/totalSampleSize+
                " 最优性分析花费"+sampleSizeOC+" 占比"+sampleSizeOC/totalSampleSize);
    }

    public double gt(int t) {

        double p =  1-alpha/((k+d-1)*1.0*(Math.log(2*t)/Math.log(2))*(Math.log(2*t)/Math.log(2)));
        TDistribution td = new TDistribution(t-1);

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
}
