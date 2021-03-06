#include "FactoredMDP.h"
#include <dai/jtree.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "XML.h"

using namespace dai;
using namespace std;


// load a factored MDP from a file
FactoredMDP::FactoredMDP(const char* filename)
{
  m_directory = filename;
  m_directory.erase(m_directory.rfind('/')+1);
  readFromXML(XMLNode::parseFile(filename, "mdp"));
}

// load a factored MPD from a file, modified by some experiments
FactoredMDP::FactoredMDP(const char* filename, const vector<Experiment> &exps,
			 const vector<int> exp_idxs)
{
  m_directory = filename;
  m_directory.erase(m_directory.rfind('/')+1);

  XMLNode mdp = XMLNode::parseFile(filename, "mdp");
  for (int i = 0; i < exps.size(); i++) {
    exps[i].populateNode(&mdp, exp_idxs[i]);
  }
  readFromXML(mdp);
}


// value iteration over a factored MDP
vector<Factor> FactoredMDP::valueIteration(vector<Factor> *Pi)
{
  VarSet states = m_S_set[0] | m_Z_set[0] | m_P_set[0];
  vector<Factor> Q(m_rewards.size(), Factor(m_A_set,0.0));
  vector<Factor> V(m_rewards.size(), Factor(states,0.0)),
    V2(m_rewards.size(), Factor(states, 0.0));
  Factor Q_total(m_A_set, 0.0);

  VarSet elim_vars = allVars(m_mdp) / m_A_set;
  for (int i = 0; i < 2; i++) elim_vars /= m_S_set[i] | m_Z_set[i] | m_P_set[i];
  FactorGraph r_mdp = variableElimination(m_mdp, elim_vars);

  Pi->resize(m_vi_horizon);
  for (int t = m_vi_horizon-1; t >= 0; t--) {
    if (m_verbose) cerr << "Value iteration " << t << endl;
    if (Pi != 0) (*Pi)[t] = Factor(states, 0.0);
    for (State s(states); s.valid(); s++) {
      FactorGraph r_mdp2 = addEvidence(r_mdp, s);
      Q_total.fill(0.0);
      for (State a(m_A_set); a.valid(); a++) {
	for (int i = 0; i < m_rewards.size(); i++) {
	  Q[i][a] = m_rewards[i][a(m_rewards[i].vars()) +
				 s(m_rewards[i].vars())];
	}
	
	Factor f = jointDistribution(addEvidence(r_mdp2, a));
	for (int i = 0; i < m_rewards.size(); i++) {
	  Q[i][a] += (f.p() * V[i].p()).sum();
	}
	for (int i = 0; i < m_rewards.size(); i++) {
	  Q_total[a] += Q[i][a];
	}
      }
      int act;
      factorMax(Q_total, &act);
      for (int i = 0; i < m_rewards.size(); i++) V2[i][s] = Q[i][act];
      if (Pi != 0) (*Pi)[t][s] = act;
    }
    for (int i = 0; i < m_rewards.size(); i++) V[i] = V2[i];
  }
  return V;
}


// compute policy value (decomposed along different rewards)
vector<Factor> FactoredMDP::policyValueDecomp(const vector<Factor> &Pi)
{
  VarSet states = m_S_set[0] | m_Z_set[0] | m_P_set[0];
  vector<Factor> V(m_rewards.size(), Factor(states,0.0)),
    V2(m_rewards.size(), Factor(states, 0.0));

  VarSet elim_vars = allVars(m_mdp) / m_A_set;
  for (int i = 0; i < 2; i++) elim_vars /= m_S_set[i] | m_Z_set[i] | m_P_set[i];
  FactorGraph r_mdp = variableElimination(m_mdp, elim_vars);

  for (int t = m_horizon-1; t >= 0; t--) {
    if (m_verbose) cerr << "Computing policy value " << t << endl;
    for (int i = 0; i < m_rewards.size(); i++) V2[i].fill(0.0);
    for (State s(states); s.valid(); s++) {
      FactorGraph r_mdp2 = addEvidence(r_mdp, s);
      State a(m_A_set, (int)Pi[max((int)Pi.size()-m_horizon+t,0)][s]);
      
      for (int i = 0; i < m_rewards.size(); i++) {
	V2[i][s] += m_rewards[i][a(m_rewards[i].vars())+s(m_rewards[i].vars())];
      }
      
      Factor f = jointDistribution(addEvidence(r_mdp2, a));
      for (int i = 0; i < m_rewards.size(); i++) {
	V2[i][s] += (f.p() * V[i].p()).sum();
      }
    }
    for (int i = 0; i < m_rewards.size(); i++) V[i] = V2[i];
  }
  return V;
}
  

// compute finite horizon value of a policy
Factor FactoredMDP::policyValue(const vector<Factor> &Pi)
{

  vector<Factor> V = policyValueDecomp(Pi);
  for (int i = 1; i < V.size(); i++) V[0] += V[i];
  return V[0];
}


// get an mpc policy for a given set of predictions
Factor FactoredMDP::mpcPartialPolicy(const map<Var, size_t> preds, int horizon)
{
  // set up DBN with transitions dependent on predictions
  FactorGraph r_mdp = addEvidence(m_mdp, preds);
  JTree jt(r_mdp, PropertySet("[verbose=0,updates=HUGIN]"));
  jt.init();
  jt.run();

  // set up small DBN (no prediction states)
  VarSet small_vars = m_A_set | m_S_set[0] | m_S_set[1] | m_Z_set[0]|m_Z_set[1];
  vector<Factor> small_factors;
  for (int i = 0; i < m_mdp.nrFactors(); i++) {
    if ((m_mdp.factor(i).vars() & small_vars) == m_mdp.factor(i).vars()) {
      small_factors.push_back(m_mdp.factor(i));
    }
  }
  FactorGraph s_mdp(small_factors);

  
  // find all factors that depend entirely on uncontrolled states
  vector<Factor*> s_factors;
  vector<vector<const Factor*> > r_factors;

  for (map<string,vector<vector<Var> > >::iterator i = m_P.begin();
       i != m_P.end(); i++){
    s_factors.push_back(&s_mdp.factor(findFactor(s_mdp, m_SZ[i->first][1])));
    r_factors.push_back(vector<const Factor*>());
    for (int j = 1; j < m_SZ[i->first].size(); j++)
      r_factors.back().push_back(&m_mdp.factor(findFactor(m_mdp,
							  m_SZ[i->first][j])));
  }
  

  // run value iteration
  VarSet states = m_S_set[0] | m_Z_set[0];
  Factor V(states, 0.0), V2(states), Pi(states, 0.0), Q(m_A_set);

  for (int t = horizon-1; t >= 0; t--) {
    // use prediction-dependent CPDs for uncontrollable variables
    for (int i = 0; i < s_factors.size(); i++) {
      if (t < r_factors[i].size()) {
	s_factors[i]->p() = jt.calcMarginal(r_factors[i][t]->vars()).p();
      }
    }

    // compute Bellman backup
    for (State s(states); s.valid(); s++) {
      FactorGraph s_mdp2 = addEvidence(s_mdp, s);
      Q.fill(0.0);
      for (State a(m_A_set); a.valid(); a++) {
	for (int i = 0; i < m_rewards.size(); i++) {
	  Q[a] += m_rewards[i][a(m_rewards[i].vars()) + s(m_rewards[i].vars())];
	}
	Factor f = jointDistribution(addEvidence(s_mdp2, a));
	Q[a] += (f.p() * V.p()).sum();
      }
      int act;
      V2[s] = factorMax(Q, &act);
      Pi[s] = act;
    }
    V = V2;
  }
  return Pi;
}


// get a complete MPC policy for all predictions
Factor FactoredMDP::mpcPolicy(int horizon)
{
  Factor Pi(m_S_set[0] | m_Z_set[0] | m_P_set[0], 0.0);
  Factor Pi2(m_S_set[0] | m_Z_set[0]);

  for (State p(m_P_set[0]); p.valid(); p++) {
    Pi2 = mpcPartialPolicy(p, horizon);
    for (State s(Pi2.vars()); s.valid(); s++) {
      Pi[p(Pi.vars()) + s(Pi.vars())] = Pi2[s];
    }
  }
  
  return Pi;
}


// return initial distribution of states in the first timestep
Factor FactoredMDP::initialDistribution()
{
  VarSet elim_vars = allVars(m_mdp) / (m_S_set[0] | m_Z_set[0] | m_P_set[0]);
  return jointDistribution(variableElimination(m_mdp, elim_vars));
}


// get a vector of doubles from either a file or direct text
vector<double> FactoredMDP::readDoubleFileOrVals(const XMLNode &node)
{
  vector<double> vals;
  double x;
  
  if (node.getAttribute("values") != 0) {
    istringstream iss(node.getAttribute("values"));
    iss >> x;
    while (iss.good()) {
      vals.push_back(x);
      iss >> x;
    }
  } else if (node.getAttribute("file") != 0) {
    ifstream fin;
    
    if (node.getAttribute("file")[0] != '/') {
      fin.open((m_directory + node.getAttribute("file")).c_str());
    } else {
      fin.open(node.getAttribute("file"));
    }
    
    if (!fin.good()) {
      cout << "Error: couldn't open file " << node.getAttribute("file") << endl;
      return vals;
    }
    
    fin >> x;
    while (fin.good()) {
      vals.push_back(x);
      fin >> x;
    }
  } else {
    cout << "Error: no file or values field" << endl;
  }
  return vals;
}


// read the mdp from an xml file
void FactoredMDP::readFromXML(const XMLNode &mdp)
{
  vector<Factor> factors;
  m_horizon = atoi(mdp.getAttribute("horizon"));
  m_vi_horizon = atoi(mdp.getAttribute("vi_horizon"));
  m_verbose = (mdp.getAttribute("verbose") != 0 &&
	       strcmp(mdp.getAttribute("verbose"), "true") == 0);
  int n_vars = 0;

  // actions
  for (int i = 0; i < mdp.nChildNode("action"); i++) {
    XMLNode a = mdp.getChildNode("action", i);
    string name(a.getAttribute("name"));
    m_A[name] = Var(n_vars++, atoi(a.getAttribute("n")));
  }
  
  // states
  for (int i = 0; i < mdp.nChildNode("state"); i++) {
    XMLNode s = mdp.getChildNode("state", i);
    string name(s.getAttribute("name"));
    m_SZ[name].push_back(Var(n_vars++, atoi(s.getAttribute("n"))));
  }

  // states, 2nd time through
  for (int i = 0; i < mdp.nChildNode("state"); i++) {
    XMLNode s = mdp.getChildNode("state", i);
    string name(s.getAttribute("name"));
    m_SZ[name].push_back(Var(n_vars++, m_SZ[name][0].states()));
  }

  // extension variables (parents of predicted variables)
  for (int i = 0; i < mdp.nChildNode("extension"); i++) {
    XMLNode p = mdp.getChildNode("extension", i);
    string name(p.getAttribute("node"));
    int horizon = atoi(p.getAttribute("horizon"));
    if (horizon <= 0) continue;

    if (m_SZ.count(name) == 0) {
      cout << "Error: no uncontrollable variable " << name << endl;
      return;
    }
    
    // add additional z variables
    for (int j = 0; j < horizon; j++) {
      m_SZ[name].push_back(Var(n_vars++, m_SZ[name][0].states()));
    }
  }
    

  // prediction variables
  for (int i = 0; i < mdp.nChildNode("prediction"); i++) {
    XMLNode p = mdp.getChildNode("prediction", i);
    string name(p.getAttribute("node"));
    int horizon = atoi(p.getAttribute("horizon"));
    if (horizon <= 0) continue;

    if (m_SZ.count(name) == 0) {
      cout << "Error: no uncontrollable variable " << name << endl;
      return;
    }
    
    // add additional z variables
    for (int j = 0; j < horizon; j++) {
      m_SZ[name].push_back(Var(n_vars++, m_SZ[name][0].states()));
    }
    
    // add P varibles
    m_P[name] = vector<vector<Var> >(horizon);
    for (int j = 0; j < horizon; j++) {
      for (int k = 0; k < horizon+1-j; k++) {
	m_P[name][j].push_back(Var(n_vars++, m_SZ[name][0].states()));
      }
    }
  }


  // rewards
  for (int i = 0; i < mdp.nChildNode("reward"); i++) {
    XMLNode r = mdp.getChildNode("reward", i);
    vector<double> vals = readDoubleFileOrVals(r);

    VarSet scope;
    for (int j = 0; j < r.nChildNode("parent"); j++) {
      if (m_A.count(r.getChildNode("parent", j).getAttribute("name")) > 0) {
	scope |= m_A[r.getChildNode("parent", j).getAttribute("name")];
      } else {
	scope |= m_SZ[r.getChildNode("parent", j).getAttribute("name")][0];
      }
    }
    m_rewards.push_back(Factor(scope, vals));
    if (r.getAttribute("scale") != 0) {
      m_reward_scales.push_back(atof(r.getAttribute("scale")));
      m_rewards.back() *= m_reward_scales.back();
    } else {
      m_reward_scales.push_back(1.0);
    }
  }


  // transition cpds
  for (int i = 0; i < mdp.nChildNode("cpd"); i++) {
    XMLNode cpd = mdp.getChildNode("cpd", i);
    vector<double> vals = readDoubleFileOrVals(cpd);

    if (stringEndsWith(cpd.getAttribute("node"), "_n")) {
      // state transition CPDs
      string name(stringStripEnd(cpd.getAttribute("node"), "_n"));
      
      if (m_SZ.count(name) == 0) {
	cout << "Error: no node named " << name << endl;
	return;
      }
      
      // populate all transition CPDs for the variable
      for (int j = 0; j < m_SZ[name].size()-1; j++) {
	VarSet scope(m_SZ[name][j+1]);
	
	for (int k = 0; k < cpd.nChildNode("parent"); k++) {
	  string pname = cpd.getChildNode("parent",k).getAttribute("name");
	  if (stringEndsWith(pname, "_n")) {
	    scope |= m_SZ[stringStripEnd(pname, "_n")][j+1];
	  } else {
	    if (m_SZ.count(pname) > 0) {
	      scope |= m_SZ[pname][j];
	    } else {
	      scope |= m_A[pname];
	    }
	  }
	}
	factors.push_back(Factor(scope, vals));
      }
    } else {
      // initial state distribution
      VarSet scope(m_SZ[cpd.getAttribute("node")][0]);
      for (int j = 0; j < cpd.nChildNode("parent"); j++) {
	scope |= m_SZ[cpd.getChildNode("parent",j).getAttribute("name")][0];
      }
      factors.push_back(Factor(scope, vals));
    }
  }

  // prediction cpds
  for (int i = 0; i < mdp.nChildNode("prediction"); i++) {
    XMLNode p = mdp.getChildNode("prediction", i);
    string name(p.getAttribute("node"));
    int horizon = atoi(p.getAttribute("horizon"));
    if (horizon <= 0) continue;
    
    for (int j = 0; j < p.nChildNode("prediction_cpd"); j++) {
      XMLNode pcpd = p.getChildNode("prediction_cpd", j);
      vector<double> vals = readDoubleFileOrVals(pcpd);

      // assign first layer of prediction variable cpds
      if (strcmp(pcpd.getAttribute("horizon"), "all") == 0 ||
	  strcmp(pcpd.getAttribute("horizon"), "1") == 0) {
	for (int k = 0; k < horizon+1; k++) {
	  factors.push_back(Factor(VarSet(m_SZ[name][k+1], m_P[name][0][k]),
				   vals));
	}
      }

      // assign subsequent layer of prediction variable cpds
      if (strcmp(pcpd.getAttribute("horizon"), "all") == 0) {
	for (int k = 1; k < horizon; k++) {
	  for (int l = 0; l < horizon-k+1; l++) {
	    factors.push_back(Factor(VarSet(m_P[name][k-1][l+1],
					    m_P[name][k][l]), vals));
	  }
	}
      } else {
	int k = atoi(pcpd.getAttribute("horizon")) - 1;
	for (int l = 0; l < horizon-k+1; l++) {
	  factors.push_back(Factor(VarSet(m_P[name][k-1][l+1],
					  m_P[name][k][l]), vals));
	}
      }
    }
  }

  
  // populate the state, action, and prediction VarSets
  for (map<string,Var>::iterator i = m_A.begin(); i != m_A.end(); i++) {
    m_A_set |= i->second;
  }

  m_S_set.resize(2);
  m_Z_set.resize(2);
  m_P_set.resize(2);
  
  for (map<string,vector<Var> >::iterator i=m_SZ.begin(); i!=m_SZ.end(); i++) {
    if (m_P.count(i->first) > 0) {
      if (i->second.size() > m_Z_set.size()) m_Z_set.resize(i->second.size());
      for (int j = 0; j < i->second.size(); j++) m_Z_set[j] |= i->second[j];

      if (m_P_set.size() < m_Z_set.size()-1) m_P_set.resize(m_Z_set.size()-1);
      for (int j = 0; j < m_P[i->first].size(); j++) {
	for (int k = 0; k < m_P[i->first][j].size(); k++) {
	  m_P_set[k] |= m_P[i->first][j][k];
	}
      }
    } else {
      if (i->second.size() > m_S_set.size()) m_S_set.resize(i->second.size());
      for (int j = 0; j < i->second.size(); j++) m_S_set[j] |= i->second[j];
    }
  }

  m_mdp = FactorGraph(factors);
}  
  


// simple variable elimination (eliminate in their default order)
FactorGraph variableElimination(const FactorGraph &fg, const VarSet &elim)
{
  vector<Factor> factors = fg.factors();
  for (VarSet::const_iterator i = elim.begin(); i != elim.end(); i++) {
    Factor f;
    for (vector<Factor>::iterator j = factors.begin(); j!=factors.end(); j++) {
      if (j->vars().contains(*i)) {
	f *= *j;
	j = factors.erase(j)-1;
      }
    }
    factors.push_back(f.marginal(f.vars() / *i));
  }
  return FactorGraph(factors);
}


// add evidence to a Factor Graph
FactorGraph addEvidence(const FactorGraph &fg, const map<Var, size_t> &evid)
{
  vector<Factor> factors = fg.factors();
  for (vector<Factor>::iterator i = factors.begin(); i != factors.end(); i++) {
    for (map<Var,size_t>::const_iterator j = evid.begin(); j!=evid.end(); j++) {
      if (i->vars().contains(j->first)) {
	*i = i->slice(VarSet(j->first), j->second);
      }
    }
  }

  return FactorGraph(factors);
}

// return joint distirbution of a factor graph
Factor jointDistribution(const FactorGraph &fg)
{
  Factor m;
  for (vector<Factor>::const_iterator i = fg.factors().begin();
       i != fg.factors().end(); i++) {
    m *= *i;
  }
  if (m.maxAbs() > 0.0) {
    m.normalize();
  }
  return m;
}



// get var set for all factors in a graph
VarSet allVars(const FactorGraph &fg)
{
  VarSet vs;

  for (vector<Factor>::const_iterator i = fg.factors().begin();
       i != fg.factors().end(); i++) vs |= i->vars();
  return vs;
}


// get maximum value and index from factor
double factorMax(const Factor &f, int *max_i)
{
  double max = f[0];
  *max_i = 0;
  
  for (int i = 1; i < f.states(); i++) {
    if (f[i] > max) {
      max = f[i];
      *max_i = i;
    }
  }
  return max;
}

// find first factor that depends on v
size_t findFactor(const FactorGraph &fg, const Var &v)
{
  size_t i;
  for (i = 0; i < fg.nrFactors(); i++) {
    if ((fg.factor(i).vars() & VarSet(v)) == VarSet(v)) {
      return i;
    }
  }
  return i;
}
  


// quick utility, see if a string end with another string
bool stringEndsWith(const string &s1, const string &s2)
{
  size_t i = s1.rfind(s2);
  return ((i != string::npos) && (i == (s1.length() - s2.length())));
}

// strip an ending from a string if it exists
string stringStripEnd(const string &s1, const string &s2)
{
  if (stringEndsWith(s1, s2)) {
    return s1.substr(0, s1.length() - s2.length());
  } else {
    return s1;
  }
}
