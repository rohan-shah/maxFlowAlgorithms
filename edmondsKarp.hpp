#ifndef EDMONDS_KARP_HEADER_GUARD
#define EDMONDS_KARP_HEADER_GUARD
#include <boost/graph/filtered_graph.hpp>
#include <boost/property_map/compose_property_map.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
namespace multistateTurnip
{
	template<class graphType, typename flowType> struct edmondsKarpPredicate
	{
	public:
		edmondsKarpPredicate()
			:residual(NULL), graph(NULL)
		{}
		edmondsKarpPredicate(flowType* residual, const graphType* graph)
			:residual(residual), graph(graph)
		{}
		typedef typename graphType::edge_descriptor edge_descriptor;
		bool operator()(const edge_descriptor& edge) const
		{
			int index = boost::get(boost::edge_index, *graph, edge);
			return residual[index] > 0;
		}
		flowType* residual;
		const graphType* graph;
	};
	template<class graphType, typename flowType> struct edmondsKarpMaxFlowScratch
	{
		typedef boost::filtered_graph<graphType, edmondsKarpPredicate<graphType, flowType> > filteredGraphType;

		typedef typename boost::property_map<graphType, boost::vertex_index_t>::const_type vertexIndexMapType;
		typedef typename boost::property_map<graphType, boost::edge_index_t>::const_type edgeIndexMapType;
		typedef boost::compose_property_map<flowType*, edgeIndexMapType> flowMapType;
		typedef boost::compose_property_map<const flowType*, edgeIndexMapType> capacityMapType;

		typedef typename boost::property_map<filteredGraphType, boost::vertex_index_t>::const_type filteredVertexIndexMapType;
		typedef typename boost::property_map<filteredGraphType, boost::edge_index_t>::type filteredEdgeIndexMapType;
		typedef boost::iterator_property_map<typename std::vector<typename filteredGraphType::vertex_descriptor>::iterator, filteredVertexIndexMapType> filteredVertexPredecessorMapType;
		typedef boost::iterator_property_map<typename std::vector<int>::iterator, filteredVertexIndexMapType> filteredDistanceMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, filteredEdgeIndexMapType> filteredColorMapType;
		typedef boost::compose_property_map<flowType*, filteredEdgeIndexMapType> filteredEdgeMapType;


		std::vector<typename filteredGraphType::vertex_descriptor> vertexPredecessor;
		std::vector<int> distance;
		std::vector<boost::default_color_type> color;
	};
	template<class graphType, typename flowType> void edmondsKarpMaxFlow(const flowType* capacity, flowType* flow, flowType* residual, const graphType& graph, int source, int sink, flowType upperBound, edmondsKarpMaxFlowScratch<graphType, flowType>& scratch, flowType& maxFlow)
	{
		std::size_t nVertices = boost::num_vertices(graph);

		scratch.vertexPredecessor.resize(nVertices);
		scratch.distance.resize(nVertices);
		scratch.color.resize(nVertices);

		typename edmondsKarpMaxFlowScratch<graphType, flowType>::edgeIndexMapType edgeIndexMap = boost::get(boost::edge_index, graph);
		typename edmondsKarpMaxFlowScratch<graphType, flowType>::flowMapType flowMap(flow, edgeIndexMap);
		typename edmondsKarpMaxFlowScratch<graphType, flowType>::flowMapType residualMap(residual, edgeIndexMap);
		typename edmondsKarpMaxFlowScratch<graphType, flowType>::capacityMapType capacityMap(capacity, edgeIndexMap);

		edmondsKarpPredicate<graphType, flowType> predicate(residual, &graph);
		while(true)
		{
			typename graphType::out_edge_iterator current, end;
			boost::tie(current, end) = boost::out_edges(source, graph);
			if(maxFlow == upperBound)
			{
				return;
			}
			else if(maxFlow > upperBound)
			{
				//This should never happen
				throw std::runtime_error("Internal error");
			}
			else
			{
				typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredGraphType filteredGraph(graph, predicate);
				typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredVertexIndexMapType filteredVertexIndexMap = boost::get(boost::vertex_index, filteredGraph);
				typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredEdgeIndexMapType filteredEdgeIndexMap = boost::get(boost::edge_index, filteredGraph);
				typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredVertexPredecessorMapType filteredVertexPredecessorMap(scratch.vertexPredecessor.begin(), filteredVertexIndexMap);
				typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredColorMapType filteredColorMap(scratch.color.begin(), filteredEdgeIndexMap);
				typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredDistanceMapType filteredDistanceMap(scratch.distance.begin(), filteredVertexIndexMap);
				typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredEdgeMapType filteredResidualMap(residual, filteredEdgeIndexMap);

				boost::dijkstra_shortest_paths(filteredGraph, source, boost::predecessor_map(filteredVertexPredecessorMap).distance_map(filteredDistanceMap).weight_map(filteredResidualMap).color_map(filteredColorMap));
				//compute capacity of bottleneck
				typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredGraphType::vertex_descriptor current = sink;
				flowType bottleneck = std::numeric_limits<flowType>::max();
				int bottleneckEdgeIndex = -1;
				if((int)scratch.vertexPredecessor[sink] == (int)sink)
				{
					return;
				}
				while((int)current != source)
				{
					typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredGraphType::vertex_descriptor next = scratch.vertexPredecessor[current];
					typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredGraphType::edge_descriptor filteredEdge = boost::edge(next, current, filteredGraph).first;
					int candidateEdgeIndex = boost::get(boost::edge_index, filteredGraph, filteredEdge);
					if(residualMap[filteredEdge] < bottleneck)
					{
						bottleneckEdgeIndex = candidateEdgeIndex;
						bottleneck = residual[bottleneckEdgeIndex];
					}
					current = next;
				}
				if(maxFlow + bottleneck > upperBound)
				{
					bottleneck = upperBound - maxFlow;
					maxFlow = upperBound;
				}
				else maxFlow += bottleneck;
				if(bottleneck < 0) throw std::runtime_error("Internal error");
				current = sink;
				while((int)current != source)
				{
					typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredGraphType::vertex_descriptor next = scratch.vertexPredecessor[current];
					typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredGraphType::edge_descriptor filteredEdge = boost::edge(next, current, filteredGraph).first;
					flowMap[filteredEdge] += bottleneck;
					if(fabs(flowMap[filteredEdge] - capacityMap[filteredEdge]) < 1e-8) flowMap[filteredEdge] = capacityMap[filteredEdge];
					residualMap[filteredEdge] -= bottleneck;

					typename edmondsKarpMaxFlowScratch<graphType, flowType>::filteredGraphType::edge_descriptor reverseEdge = boost::get(boost::edge_reverse, filteredGraph, filteredEdge);
					flowMap[reverseEdge] = -flowMap[filteredEdge];
					residualMap[reverseEdge] = capacityMap[reverseEdge] - flowMap[reverseEdge];
					current = scratch.vertexPredecessor[current];
				}
			}
		}
	}
}
#endif

