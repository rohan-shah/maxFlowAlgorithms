#ifndef UPDATE_MAX_FLOW_INCREMENTAL_HEADER_GUARD
#define UPDATE_MAX_FLOW_INCREMENTAL_HEADER_GUARD
#include "edmondsKarp.hpp"
namespace multistateTurnip
{
	template<class graphType> struct updateMaxFlowIncrementalArgs
	{
		updateMaxFlowIncrementalArgs(const graphType& graph, edmondsKarpMaxFlowScratch& scratch)
			:graph(graph), scratch(scratch)
		{}
		double* capacity, *flow, *residual;
		graphType::edge_descriptor edge;
		const graphType& graph;
		double newCapacity;
		edmondsKarpMaxFlowScratch& scratch;
		int source, sink;
		std::vector<double> working1, working2, working3;
		int nDirectedEdges;
	};
	template<class graphType> void updateMaxFlowIncremental(updateMaxFlowIncrementalArgs<graphType>& args, double previousMaxFlow, double& newMaxFlow)
	{
		int edgeIndex = boost::get(boost::edge_index, args.graph, args.edge);
		graphType::edge_descriptor reverseEdge = boost::get(boost::edge_reverse, args.graph, args.edge);
		int reverseEdgeIndex = boost::get(boost::edge_index, args.graph, reverseEdge);
		if(args.flow[edgeIndex] < 0)
		{
			std::swap(edgeIndex, reverseEdgeIndex);
			std::swap(reverseEdge, args.edge);
		}
		if(args.newCapacity > args.capacity[edgeIndex])
		{
			if(args.flow[edgeIndex] < args.capacity[edgeIndex])
			{
				newMaxFlow = previousMaxFlow;
				return;
			}
			else
			{
				args.working1.resize(args.nDirectedEdges);
				args.working2.resize(args.nDirectedEdges);
				args.working3.resize(args.nDirectedEdges);
				memcpy(&(args.working1.front()), args.capacity, sizeof(double)*args.nDirectedEdges);
				memcpy(&(args.working2.front()), args.flow, sizeof(double)*args.nDirectedEdges);
				memcpy(&(args.working3.front()), args.residual, sizeof(double)*args.nDirectedEdges);
				//increase capacity
				args.working1[edgeIndex] = args.working1[reverseEdgeIndex] = args.newCapacity;
				//update residuals
				args.working3[edgeIndex] = args.working1[edgeIndex] - args.working2[edgeIndex];
				args.working3[reverseEdgeIndex] = args.working1[reverseEdgeIndex] - args.working2[reverseEdgeIndex];
				//update max flow
				newMaxFlow = previousMaxFlow;
				edmondsKarpMaxFlow(&(args.working1[0]), &(args.working2[0]), &(args.working3[0]), args.graph, args.source, args.sink, std::numeric_limits<double>::infinity(), args.scratch, newMaxFlow);
			}
		}
		else
		{
			if(args.flow[edgeIndex] <= args.newCapacity)
			{
				newMaxFlow = previousMaxFlow;
				return;
			}
			else
			{
				double excess = args.flow[edgeIndex] - args.newCapacity;
				args.working1.resize(args.nDirectedEdges);
				args.working2.resize(args.nDirectedEdges);
				args.working3.resize(args.nDirectedEdges);
				//working1 is the new residual matrix
				memcpy(&(args.working1.front()), args.residual, sizeof(double)*args.nDirectedEdges);
				args.working1[edgeIndex] = 0;
				args.working1[reverseEdgeIndex] = 0;
				//working2 is a copy of the new residual matrix
				memcpy(&(args.working2.front()), &(args.working1.front()), sizeof(double)*args.nDirectedEdges);
				//Outputs
				std::fill(args.working3.begin(), args.working3.end(), 0);
				double rerouted = 0;
				edmondsKarpMaxFlow(&(args.working1[0]), &(args.working3[0]), &(args.working2[0]), args.graph, args.edge.m_source, args.edge.m_target, excess, args.scratch, rerouted);
				newMaxFlow = previousMaxFlow - excess + rerouted;
			}
		}
	}
}
#endif
