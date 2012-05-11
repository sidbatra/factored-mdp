#include "Experiment.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cmath>

using namespace std;

// constructor just load xml file
Experiment::Experiment(const char *filename)
{
  m_experiment = XMLNode::parseFile(filename);
}

// output the number of different runs
int Experiment::numSteps() const
{
  if (m_experiment.getAttribute("steps") != NULL) {
    return atoi(m_experiment.getAttribute("steps"));
  } else {
    if (m_experiment.getAttribute("interp") == 0 ||
	strcmp(m_experiment.getAttribute("interp"), "linear") == 0) {
      return floor((atof(m_experiment.getAttribute("max")) -
		    atof(m_experiment.getAttribute("min"))) /
		   atof(m_experiment.getAttribute("step")) + 1e-10) + 1;
    } else if (strcmp(m_experiment.getAttribute("interp"), "exponential") == 0) {
      return floor(log(atof(m_experiment.getAttribute("max")) /
		       atof(m_experiment.getAttribute("min"))) /
		   log(atof(m_experiment.getAttribute("step"))) + 1e-10) + 1;
    }
  }
}


// alter an xml tree to change items corresponding to the experiment
bool Experiment::populateNode(XMLNode *root, int trial_idx) const
{
  // search until we find the node of interest
  XMLNode e_node = m_experiment.getChildNode();
  XMLNode node = *root;
  if (strcmp(e_node.getName(), node.getName()) != 0) return false;
    
  while (e_node.nChildNode() > 0) {
    e_node = e_node.getChildNode();
    if (e_node.nAttribute() > 0) {
      node = node.getChildNodeWithAttribute(e_node.getName(),
					    e_node.getAttributeName(),
					    e_node.getAttributeValue());
    } else {
      node = node.getChildNode(e_node.getName());
    }
  }

  // compute linear interpolate amounto
  double alpha;

  if (numSteps() > 1) {
    if (m_experiment.getAttribute("interp") == 0 ||
	strcmp(m_experiment.getAttribute("interp"), "linear") == 0) {
      alpha = (double)trial_idx / (numSteps() - 1);
    } else {
      // convert exponential fraction to linear fraction
      double r = exp(log(atof(m_experiment.getAttribute("max")) /
			 atof(m_experiment.getAttribute("min"))) / (numSteps()-1));
      double val = atof(m_experiment.getAttribute("min")) * pow(r, trial_idx);
      if (val > atof(m_experiment.getAttribute("min")))
	val = atof(m_experiment.getAttribute("min"));
      alpha = (val - atof(m_experiment.getAttribute("min"))) /
	(atof(m_experiment.getAttribute("max")) -
	 atof(m_experiment.getAttribute("min")));
    }
  } else {
    alpha = 0.0;
  }
  
  cout << (1-alpha)*atof(m_experiment.getAttribute("min")) +
    alpha * atof(m_experiment.getAttribute("max")) << " ";
  
  if (strcmp(m_experiment.getAttribute("type"), "integer") == 0) {
    node.updateAttribute(getIntegerString(alpha).c_str(), NULL,
			 m_experiment.getAttribute("attribute"));
  } else if (strcmp(m_experiment.getAttribute("type"), "double") == 0) {
    node.updateAttribute(getDoubleString(alpha).c_str(), NULL,
			 m_experiment.getAttribute("attribute"));
  } else if (strcmp(m_experiment.getAttribute("type"), "vector") == 0) {
    node.updateAttribute(getVectorString(alpha).c_str(), NULL,
			 m_experiment.getAttribute("attribute"));
  } else if (strcmp(m_experiment.getAttribute("type"), "file") == 0) {
    node.updateAttribute(getFileString(alpha).c_str(), NULL,
			 m_experiment.getAttribute("attribute"));
  } else {
    cout << "Error: not a valid type "
	 << m_experiment.getAttribute("type") << endl;
    return false;
  }

  root->writeToFile("xml.out");

  return true;
}


// interpolate integers
string Experiment::getIntegerString(double alpha) const
{
  ostringstream oss;
  oss << (int)((1-alpha)*atoi(m_experiment.getAttribute("min")) +
	       alpha*atoi(m_experiment.getAttribute("max")) + 1e-10);
  return string(oss.str());
}

// interpolate doubles
string Experiment::getDoubleString(double alpha) const
{
  ostringstream oss;
  oss << (double)((1-alpha)*atof(m_experiment.getAttribute("min")) +
		  alpha*atof(m_experiment.getAttribute("max")));
  return string(oss.str());
}

// interpolate vector of doubles
string Experiment::getVectorString(double alpha) const
{
  istringstream iss_min(m_experiment.getAttribute("min"));
  istringstream iss_max(m_experiment.getAttribute("max"));
  ostringstream oss;
  double x_min, x_max;
  
  iss_min >> x_min; iss_max >> x_max;
  while (iss_min.good()) {
    oss << (1-alpha)*x_min + alpha*x_max << " ";
    iss_min >> x_min; iss_max >> x_max;
  }

  return oss.str();
}
    

// interpolate files (create a temp file with interpolated values)
string Experiment::getFileString(double alpha) const
{
  ifstream fin_min(m_experiment.getAttribute("min"));
  ifstream fin_max(m_experiment.getAttribute("max"));
  double x_min, x_max;
  
  char filename[256] = "/tmp/interp_file_XXXXXX";
  int fd = mkstemp(filename);
  ofstream fout(filename);

  close(fd);
  fin_min >> x_min; fin_max >> x_max;
  while (fin_min.good()) {
    fout << (1-alpha)*x_min + alpha*x_max << endl;
    fin_min >> x_min; fin_max >> x_max;
  }
  
  return string(filename);
}
