#ifndef UPDATE_FLOW_INCREMENTAL_HEADER_GUARD
#define UPDATE_FLOW_INCREMENTAL_HEADER_GUARD
#include "edmondsKarp.hpp"
namespace multistateTurnip
{
	template<class graphType, typename flowType> struct updateFlowIncrementalArgs
	{
		updateFlowIncrementalArgs(const graphType& graph, edmondsKarpMaxFlowScratch<graphType, flowType>& scratch)
			:graph(graph), scratch(scratch)
		{}
		flowType* capacity, *flow, *residual;
		typename graphType::edge_descriptor edge;
		const graphType& graph;
		flowType newCapacity;
		edmondsKarpMaxFlowScratch<graphType, flowType>& scratch;
		int source, sink;
		std::vector<flowType> working1, working2, working3;
		int nDirectedEdges;
	};
	template<class graphType, typename flowType> void updateFlowIncremental(updateFlowIncrementalArgs<graphType, flowType>& args, flowType previousMaxFlow, flowType& newMaxFlow)
	{
		int edgeIndex = boost::get(boost::edge_index, args.graph, args.edge);
		typename graphType::edge_descriptor reverseEdge = boost::get(boost::edge_reverse, args.graph, args.edge);
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
				args.capacity[edgeIndex] = args.capacity[reverseEdgeIndex] = args.newCapacity;
				args.residual[edgeIndex] = args.newCapacity - args.flow[edgeIndex];
				args.residual[reverseEdgeIndex] = args.newCapacity - args.flow[reverseEdgeIndex];
				newMaxFlow = previousMaxFlow;
				return;
			}
			else
			{
				args.capacity[edgeIndex] = args.capacity[reverseEdgeIndex] = args.newCapacity;
				args.residual[edgeIndex] = args.newCapacity - args.flow[edgeIndex];
				args.residual[reverseEdgeIndex] = args.newCapacity - args.flow[reverseEdgeIndex];
				newMaxFlow = previousMaxFlow;
				edmondsKarpMaxFlow<graphType, flowType>(args.capacity, args.flow, args.residual, args.graph, args.source, args.sink, std::numeric_limits<flowType>::max(), args.scratch, newMaxFlow);
			}
		}
		else
		{
			if(args.flow[edgeIndex] <= args.newCapacity)
			{
				newMaxFlow = previousMaxFlow;
				args.capacity[edgeIndex] = args.capacity[reverseEdgeIndex] = args.newCapacity;
				args.residual[edgeIndex] = args.newCapacity - args.flow[edgeIndex];
				args.residual[reverseEdgeIndex] = args.newCapacity - args.flow[reverseEdgeIndex];
				return;
			}
			else
			{
				flowType excess = args.flow[edgeIndex] - args.newCapacity;
				args.capacity[edgeIndex] = args.capacity[reverseEdgeIndex] = args.newCapacity;
				args.working1.resize(args.nDirectedEdges);
				args.working2.resize(args.nDirectedEdges);
				args.working3.resize(args.nDirectedEdges);

				//working1 is the new residual matrix, after we alter the flow so that the edge is saturated
				memcpy(&(args.working1.front()), args.residual, sizeof(flowType)*args.nDirectedEdges);
				args.working1[edgeIndex] = 0;
				args.working1[reverseEdgeIndex] = 0;//2 * args.newCapacity;
				//working2 is a copy of the new residual matrix
				memcpy(&(args.working2.front()), &(args.working1.front()), sizeof(flowType)*args.nDirectedEdges);
				//Outputs
				std::fill(args.working3.begin(), args.working3.end(), 0);
				flowType rerouted = 0;
				edmondsKarpMaxFlow<graphType, flowType>(&(args.working1[0]), &(args.working3[0]), &(args.working2[0]), args.graph, args.edge.m_source, args.edge.m_target, excess, args.scratch, rerouted);
				//Update flow and residual
				for(int i = 0; i < args.nDirectedEdges; i++) 
				{
					args.flow[i] += args.working3[i];
					args.residual[i] -= args.working3[i];
				}
				args.residual[edgeIndex] = 0;
				args.residual[reverseEdgeIndex] = 0;
				if((int)args.edge.m_source != args.source)
				{
					memcpy(&(args.working2[0]), &(args.residual[0]), sizeof(flowType)*args.nDirectedEdges);
					memcpy(&(args.working1[0]), &(args.residual[0]), sizeof(flowType)*args.nDirectedEdges);
					std::fill(args.working3.begin(), args.working3.end(), 0);
					flowType unused = 0;
					edmondsKarpMaxFlow(&(args.working1[0]), &(args.working3[0]), &(args.working2[0]), args.graph, args.edge.m_source, args.source, excess - rerouted, args.scratch, unused);
					for(int i = 0; i < args.nDirectedEdges; i++) args.flow[i] += args.working3[i];
				}
				if((int)args.edge.m_target != args.sink)
				{
					memcpy(&(args.working2[0]), &(args.residual[0]), sizeof(flowType)*args.nDirectedEdges);
					memcpy(&(args.working1[0]), &(args.residual[0]), sizeof(flowType)*args.nDirectedEdges);
					std::fill(args.working3.begin(), args.working3.end(), 0);
					flowType unused = 0;
					edmondsKarpMaxFlow(&(args.working1[0]), &(args.working3[0]), &(args.working2[0]), args.graph, args.sink, args.edge.m_target, excess - rerouted, args.scratch, unused);
					for(int i = 0; i < args.nDirectedEdges; i++) args.flow[i] += args.working3[i];
				}
				newMaxFlow = previousMaxFlow - excess + rerouted;
				//update residual
				args.flow[edgeIndex] = args.newCapacity;
				args.flow[reverseEdgeIndex] = -args.newCapacity;
				for(int i = 0; i < args.nDirectedEdges; i++) args.residual[i] = args.capacity[i] - args.flow[i];
			}
		}
	}

}
#endif
