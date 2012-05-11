#ifndef FACTORED_MDP_H
#define FACTORED_MDP_H

/*
  A collection of routines for loading / solving factored MDPs using
  the libdai library.  We're also including some simple "efficient"
  inference and factor graph manipulation routines, because the libdai
  built-in exact inference methods are _slow_.
*/

#include <dai/factor.h>
#include <dai/factorgraph.h>
#include <dai/var.h>
#include <dai/varset.h>
#include <vector>
#include <map>
#include <string>
#include "XML.h"
#include "Experiment.h"

using namespace dai;


  class FactoredMDP {
  public:
    FactoredMDP(const char* filename);
    FactoredMDP(const char* filename, const std::vector<Experiment> &exps,
		const std::vector<int> exp_idxs);
    ~FactoredMDP() {};

    std::vector<Factor> valueIteration(std::vector<Factor> *Pi = 0);
    Factor policyValue(const std::vector<Factor> &Pi);
    std::vector<Factor> policyValueDecomp(const std::vector<Factor> &Pi);
    
    Factor mpcPartialPolicy(const std::map<Var, size_t> preds, int horizon);
    Factor mpcPolicy(int horizon);
    Factor initialDistribution();

    int horizon() { return m_horizon; }
    int viHorizon() { return m_vi_horizon; }
    double rewardScale(int i) { return m_reward_scales[i]; }

    Var stateVar(std::string name) { return m_SZ[name][0]; }
    Var predVar(std::string name, int h) { return m_P[name][h][0]; }
    VarSet allStates() { return m_S_set[0] | m_Z_set[0] | m_P_set[0]; }

  private:
    std::vector<double> readDoubleFileOrVals(const XMLNode &node);
    void readFromXML(const XMLNode &mdp);
    
    FactorGraph m_mdp;
    std::vector<Factor> m_rewards;
    std::vector<double> m_reward_scales;
    int m_horizon, m_vi_horizon, m_n_vars;
    std::string m_directory;
    bool m_verbose;

    std::map<std::string, Var> m_A;
    std::map<std::string, std::vector<Var> > m_SZ;
    std::map<std::string, std::vector<std::vector<Var> > > m_P;

    VarSet m_A_set;
    std::vector<VarSet> m_S_set;
    std::vector<VarSet> m_Z_set;
    std::vector<VarSet> m_P_set;
  };



  // general factor graph utility functions  
  FactorGraph variableElimination(const FactorGraph &fg, const VarSet &elim);
  FactorGraph addEvidence(const FactorGraph &fg,
			  const std::map<Var, size_t> &evid);
  Factor jointDistribution(const FactorGraph &fg);
  VarSet allVars(const FactorGraph &fg);
  double factorMax(const Factor &f, int *max_i);
  size_t findFactor(const FactorGraph &fg, const Var &v);

  // simple string utils for reading xml file
  bool stringEndsWith(const std::string &s1, const std::string &s2);
  std::string stringStripEnd(const std::string &s1, const std::string &s2);


#endif

