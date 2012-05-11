#include "FactoredMDP.h"
#include "Experiment.h"
#include <dai/alldai.h>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstdio>

using namespace std;
using namespace dai;


int main(int argc, char* argv[])
{
  if (argc < 2) {
    cout << "Usage: ./main <mdp.xml> <num_experiments> <exp_1.xml> ..." << endl;
    return 1;
  }

  vector<Experiment> exps;
  vector<size_t> sizes;

  if (argc > 2 && atoi(argv[2]) > 0) {
    for (int i = 0; i < atoi(argv[2]); i++) {
      exps.push_back(Experiment(argv[3+i]));
      sizes.push_back(exps.back().numSteps());
    } 
  }

  int idx = 0;
  for (multifor e(sizes); e.valid(); e++) {
    vector<int> e_idxs;
    FactoredMDP *mdp;
    
    if (argc > 2 && atoi(argv[2]) > 0) {
      for (int i = 0; i < atoi(argv[2]); i++) e_idxs.push_back(e[i]);
      mdp = new FactoredMDP(argv[1], exps, e_idxs);
      //for (int i = 0; i < atoi(argv[2]); i++) cout << e[i] << " ";
    } else {
      mdp = new FactoredMDP(argv[1]);
    }

    // compute value and print the expected value
    Factor d0 = mdp->initialDistribution();
    vector<Factor> Pi;
    vector<Factor> V = mdp->valueIteration(&Pi);
    for (int i = 0; i < V.size(); i++) {
      cout << setprecision(10)
	   << (d0 * V[i]).sum()/mdp->rewardScale(i)/mdp->viHorizon() << " ";
    }
    cout << endl;
    
    // run mpc
    /*
    vector<Factor> Pi_mpc(mdp->viHorizon());
    for (int i = 0; i < mdp->viHorizon(); i++) {
      Pi_mpc[i] = mdp->mpcPolicy(mdp->viHorizon()-i);
    }
    vector<Factor> V_mpc = mdp->policyValueDecomp(Pi_mpc);
    for (int i = 0; i < V.size(); i++) {
      cout << setprecision(10) << (d0 * V_mpc[i]).sum() << " ";
    }
    cout << endl;
    */
    
    // output policy
    char filename[256];
    sprintf(filename, "policies/Pi_%d.txt", idx++);
    ofstream fout(filename);
    for (int i = 0; i < Pi[0].states(); i++) {
      fout << Pi[0][i] << endl;
    }
    fout.close();

    delete mdp;
  }

  return 0;
}

