#pragma once

#include <vector>
using namespace std;

#include "gurobi_c++.h"

#include "instance.hpp"
#include "common.hpp"
#include "solution.hpp"

class CutModel
{
private:
    Instance &I;
    Solution &S;
    GRBEnv env;
    GRBModel model;
    map<VERTEXID, GRBVar> x; // Cut set variables
    map<VERTEXID, GRBVar> y; // Reachability variables

public:
    CutModel(Solution& _S, TIME threshold, CAPACITY num_resources, unsigned timelimit=300) : I(_S.I), S(_S), env(), model(env)
    {
        // Initialize model
        model.set(GRB_StringAttr_ModelName, "MinimizeReachableVertices");
        model.set(GRB_IntParam_OutputFlag, true);
        model.set(GRB_DoubleParam_TimeLimit, timelimit);
        model.set(GRB_IntParam_Threads, 1);
        // Find source vertices
        set<VERTEXID> sources;
        for (VERTEXID v : I.get_V())
            if (S.fire_arrival_time[v] < threshold)
                sources.insert(v);
        // Create binary variables x_i (cut indicator) and y_i (reachability indicator)
        for (VERTEXID v : I.get_V())
        {
            x[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + std::to_string(v));
            y[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + std::to_string(v));
        }
        // Objective: Minimize sum of y_i (reachable vertices)
        GRBLinExpr objective = 0;
        for (VERTEXID v : I.get_V())
        {
            objective += y[v];
        }
        // Objective function
        model.setObjective(objective, GRB_MINIMIZE);
        //  Limit the cut size to num_resources
        GRBLinExpr cut_size_constraint = 0;
        for (VERTEXID v : I.get_V())
        {
            cut_size_constraint += x[v];
        }
        model.addConstr(cut_size_constraint <= num_resources, "CutSizeConstraint");
        // Source vertices in S must be reachable
        for (VERTEXID s : sources)
        {
            model.addConstr(y[s] == 1, "SourceReachability_" + std::to_string(s));
        }
        // Reachability through arcs
        for (VERTEXID u : I.get_V())
        {
            for (auto &[v, t_uv] : I.get_outgoing_arcs(u))
            {
                if (y.find(v) != y.end())
                    model.addConstr(y[u] <= y[v] + x[u], "ReachabilityPropagation_" + std::to_string(u) + "_" + std::to_string(v));
            }
        }
        // Ensure source vertices cannot be in the cut set
        for (VERTEXID s : sources)
        {
            model.addConstr(x[s] == 0, "SourceNotInCut_" + std::to_string(s));
        }
        model.update();
    }


    bool optimize(vector<VERTEXID>& P)
    {
        model.optimize();
        P.clear();
        if (model.get(GRB_IntAttr_SolCount) > 0)
        {
            for(auto& [v, var] : x)
                if (x[v].get(GRB_DoubleAttr_X) > 0.5)
                    P.push_back(v);
            return true;
        }
        else
        {
            return false;
        }
    }
};
