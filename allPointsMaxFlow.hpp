#ifndef ALL_POINTS_MAX_FLOW_HEADER_GUARD
#define ALL_POINTS_MAX_FLOW_HEADER_GUARD
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include "edmondsKarp.hpp"
#include <boost/graph/undirected_dfs.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include "custom_dfs.hpp"
namespace allPointsMaxFlow
{
	//Predicate to filter out edges with zero residual flow
	template<class inputGraph, typename flowType = double> struct flowPredicate
	{
	public:
		typedef typename inputGraph::edge_descriptor edge_descriptor;
		flowPredicate()
			:graph(NULL)
		{}
		flowPredicate(const inputGraph& graph, typename std::vector<flowType>::const_iterator capacities, typename std::vector<flowType>::const_iterator residualCapacities)
			: graph(&graph), capacities(capacities), residualCapacities(residualCapacities)
		{}
		bool operator()(const edge_descriptor& edge) const
		{
			int index = boost::get(boost::edge_index, *graph, edge);
			return *(residualCapacities + index) > 0;// || (*residualCapacities)[reverseIndex] < (*capacities)[reverseIndex];
		}
	private:
		const inputGraph* graph;
		typename std::vector<flowType>::const_iterator capacities;
		typename std::vector<flowType>::const_iterator residualCapacities;
	};
	//Working data for the all points max flow code. Includes loads of typedefs
	template<class inputGraph, typename flowType = double> struct allPointsMaxFlowScratch
	{
		//Working data for the edmonds Karp call. 
		std::vector<flowType> edgeResidualCapacityVector;
		//These two color vectors are used in depth-first searches of the flow-equivalent tree (to extract the flows) and of the flow-equivalent tree (in the process of actually constructing it)
		std::vector<boost::default_color_type> colorVector1;
		std::vector<boost::default_color_type> colorVector2;

		//The type of the flow equivalent tree. I'm calling this a star graph because during construction of the flow-equivalent tree, this graph starts out as a star, and remains a tree). 
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_index_t, int, boost::property<boost::edge_weight_t, flowType> > > starGraphType;
		//Typedefs for depth first searches of the flow-equivalent tree
		typedef typename boost::property_map<starGraphType, boost::edge_index_t>::const_type starEdgeIndexMapType;
		typedef typename boost::property_map<starGraphType, boost::vertex_index_t>::const_type starVertexIndexMapType;
		typedef typename boost::property_map<starGraphType, boost::edge_weight_t>::const_type starEdgeWeightMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, starVertexIndexMapType> starVertexColorMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, starEdgeIndexMapType> starEdgeColorMapType;
		
		//Typedefs for depth-first searches of the filtered graph
		typedef boost::filtered_graph<inputGraph, flowPredicate<inputGraph, flowType> > filteredGraph;
		typedef typename boost::property_map<filteredGraph, boost::vertex_index_t>::const_type filteredVertexIndexMapType;
		typedef typename boost::property_map<filteredGraph, boost::edge_index_t>::type filteredEdgeIndexMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, filteredVertexIndexMapType> filteredVertexColorMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, filteredEdgeIndexMapType> filteredEdgeColorMapType;
		//Scratch data for the edmonds Karp calls. 
		typename multistateTurnip::edmondsKarpMaxFlowScratch<inputGraph, flowType> edmondsKarpScratch;
	};
	//This is the state for the visitor used to traverse the flow-equivalent tree and get out the all-points max flows
	template<typename inputGraph, typename flowType = double> struct minimumEdgeWeightVisitorState
	{
	public:
		minimumEdgeWeightVisitorState(typename std::vector<flowType>::iterator flowMatrix)
			: flowMatrix(flowMatrix)
		{}
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::vertex_descriptor start;
		std::vector<typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_descriptor> edges;
		std::vector<flowType> minimumValues;
		typename std::vector<flowType>::iterator flowMatrix;
		std::size_t nVertices;
	};
	//This visitor traverses the flow-equivalent ntree and gets out the all-points max flow. It does this by keeping track of the minimum capacity edge. 
	template<class inputGraph, typename flowType> class minimumEdgeWeightVisitor : public boost::default_dfs_visitor
	{
	public:
		minimumEdgeWeightVisitor(minimumEdgeWeightVisitorState<inputGraph, flowType>& stateRef)
			:stateRef(stateRef)
		{}
		void tree_edge(const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_descriptor& e, const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType& star)
		{
			flowType currentEdgeWeight = boost::get(boost::edge_weight, star, e);
			//If there are no edges at the moment, this edge must be the bottleneck
			if(stateRef.minimumValues.size() == 0)
			{
				stateRef.minimumValues.push_back(currentEdgeWeight);
				stateRef.edges.push_back(e);
			}
			//Otherwise check if this edge has lower capacity than the previous lowest capacity
			else if (currentEdgeWeight < *stateRef.minimumValues.rbegin())
			{
				stateRef.minimumValues.push_back(currentEdgeWeight);
				stateRef.edges.push_back(e);
			}
			//If we've somehow looped around and reached the starting vertex again, throw an exception. This should be impossible because the graph this visitor is applied to is a tree.
			if(e.m_target == stateRef.start) throw std::runtime_error("");
			//Put in the current bottleneck as the flow between the start vertex and the end-point of the current edge. 
			*(stateRef.flowMatrix + e.m_target + stateRef.nVertices * stateRef.start) = *(stateRef.flowMatrix + stateRef.start + stateRef.nVertices * e.m_target) = stateRef.minimumValues.back();
		}
		void finish_edge(const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_descriptor e, const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType& star)
		{
			if(stateRef.edges.size() == 0) return;
			typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_descriptor top = stateRef.edges.back();
			//If this edge is the current bottleneck, remove it from the vectors.
			if(e == top || (e.m_target == top.m_source && e.m_source == top.m_target))
			{
				stateRef.edges.pop_back();
				stateRef.minimumValues.pop_back();
			}
		}
	private:
		minimumEdgeWeightVisitorState<inputGraph, flowType>& stateRef;
	};
	template<class inputGraph, typename flowType> void extractFlow(const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType& tree, typename std::vector<flowType>::iterator& flowMatrix, allPointsMaxFlowScratch<inputGraph, flowType>& scratch)
	{
		const std::size_t nVertices = boost::num_vertices(tree);
		const std::size_t nEdges = boost::num_edges(tree);
		scratch.colorVector1.resize(nVertices);
		scratch.colorVector2.resize(nEdges);
		//If there are only two vertices in the star graph, then there is only one edge and its weight is the flow between the two vertices
		if (nVertices == 2)
		{
			typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_iterator begin, end;
			boost::tie(begin, end) = boost::edges(tree);
			*flowMatrix = *(flowMatrix + 1) = boost::get(boost::edge_weight, tree, *begin);
		}
		else if (nVertices < 2)
		{
			return;
		}
		//Property maps for the depth first search
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starVertexIndexMapType vertexIndexMap = boost::get(boost::vertex_index, tree);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starEdgeIndexMapType edgeIndexMap = boost::get(boost::edge_index, tree);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starVertexColorMapType vertexColorMap(scratch.colorVector1.begin(), vertexIndexMap);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starEdgeColorMapType edgeColorMap(scratch.colorVector2.begin(), edgeIndexMap);

		//Underlying state object for the depth first search
		minimumEdgeWeightVisitorState<typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType, flowType> visitorState(flowMatrix);
		visitorState.nVertices = nVertices;
		
		//Visitor for the depth first search
		minimumEdgeWeightVisitor<typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType, flowType> visitor(visitorState);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::vertex_iterator current, end;
		boost::tie(current, end) = boost::vertices(tree);
		//Use every vertex as a starting point
		for(; current != end; current++)
		{
			//If we didn't properly reverse every edge from the internal vectors than throw an exception. This actually gets triggered on some boost versions because there is a bug which means that finish_edge is not properly called. 
			if(visitorState.edges.size() > 0 || visitorState.minimumValues.size() > 0) throw std::runtime_error("Internal error");
			//Run a custom depth first search that calls finish_edge as it backs out along every edge
			visitorState.start = *current;
			std::fill(scratch.colorVector1.begin(), scratch.colorVector1.end(), boost::white_color);
			std::fill(scratch.colorVector2.begin(), scratch.colorVector2.end(), boost::white_color);
			boost::custom_undirected_dfs(tree, visitor, vertexColorMap, edgeColorMap, *current);
		}
	}
	template<class inputGraph, typename flowType> void allPointsMaxFlow(typename std::vector<flowType>::iterator flowMatrix, const typename std::vector<flowType>::const_iterator capacities, const inputGraph& graph, allPointsMaxFlowScratch<inputGraph, flowType>& scratch)
	{
		BOOST_STATIC_ASSERT(
			boost::is_base_and_derived<
				boost::directed_tag, 
				typename boost::graph_traits<inputGraph>::directed_category
			>::value
		); 
		const std::size_t nVertices = boost::num_vertices(graph);
		const std::size_t nEdges = boost::num_edges(graph);
		//Resize scratch data
		scratch.edgeResidualCapacityVector.resize(nEdges);
		scratch.colorVector1.resize(nVertices);
		scratch.colorVector2.resize(nEdges);

		//construct star graph
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType star(nVertices);
		for (std::size_t i = 1; i < nVertices; i++)
		{
			boost::add_edge(i, 0, (int)i-1, star);
		}
		for (std::size_t s = 1; s < nVertices; s++)
		{
			//At this point there is a single edge incident to the vertex
			typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_descriptor uniqueEdge = *boost::out_edges(s, star).first;
			typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::vertex_descriptor t = uniqueEdge.m_target;

			//Call Edgmonds Karp
			std::copy(capacities, capacities + nEdges, scratch.edgeResidualCapacityVector.begin());
			std::fill(flowMatrix, flowMatrix + nEdges, 0);
			flowType maxFlow = 0;
			multistateTurnip::edmondsKarpMaxFlow<inputGraph, flowType>(&*capacities, &*flowMatrix, &scratch.edgeResidualCapacityVector.front(), graph, s, t, std::numeric_limits<flowType>::max(), scratch.edmondsKarpScratch, maxFlow);
			boost::put(boost::edge_weight, star, uniqueEdge, maxFlow);

			//get out the two components of graph which are seperated by the mincut we just found
			typename allPointsMaxFlowScratch<inputGraph, flowType>::filteredGraph filtered(graph, flowPredicate<inputGraph, flowType>(graph, capacities, scratch.edgeResidualCapacityVector.begin()));
			typename allPointsMaxFlowScratch<inputGraph, flowType>::filteredVertexIndexMapType filteredVertexMap = boost::get(boost::vertex_index, filtered);
			typename allPointsMaxFlowScratch<inputGraph, flowType>::filteredEdgeIndexMapType filteredEdgeMap = boost::get(boost::edge_index, filtered);
			typename allPointsMaxFlowScratch<inputGraph, flowType>::filteredVertexColorMapType filteredVertexColorMap(scratch.colorVector1.begin(), filteredVertexMap);
			typename allPointsMaxFlowScratch<inputGraph, flowType>::filteredEdgeColorMapType filteredEdgeColorMap(scratch.colorVector2.begin(), filteredEdgeMap);
			std::fill(scratch.colorVector1.begin(), scratch.colorVector1.end(), boost::white_color);
			std::fill(scratch.colorVector2.begin(), scratch.colorVector2.end(), boost::white_color);
			boost::dfs_visitor<boost::null_visitor> visitor = boost::make_dfs_visitor(boost::null_visitor());
			boost::detail::undir_dfv_impl(filtered, s, visitor, filteredVertexColorMap, filteredEdgeColorMap);
			//at this point one component is white in the colour vector, and one is black. 
			for (std::size_t i = s + 1; i < nVertices; i++)
			{
				bool isNeighbour = false;
				typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::out_edge_iterator current, end;
				boost::tie(current, end) = boost::out_edges(i, star);
				while (current != end)
				{
					if (current->m_target == t) 
					{
						isNeighbour = true;
						break;
					}
					current++;
				}
				if (scratch.colorVector1[i] == boost::black_color && isNeighbour)
				{
					int oldEdgeIndex = boost::get(boost::edge_index, star, *current);
					boost::remove_edge(*current, star);
					boost::add_edge(i, s, oldEdgeIndex, star);
				}
			}
		}
		//Extract the max flow from the flow-equivalent tree (the star graph)
		std::fill(flowMatrix, flowMatrix + nEdges, std::numeric_limits<flowType>::max());
		extractFlow<inputGraph, flowType>(star, flowMatrix, scratch);
	}
	template<class inputGraph, typename flowType> void allPointsMaxFlow(std::vector<flowType>& flowMatrix, const std::vector<flowType>& capacities, const inputGraph& graph, allPointsMaxFlowScratch<inputGraph, flowType>& scratch)
	{
		allPointsMaxFlow<inputGraph, flowType>(flowMatrix.begin(), capacities.begin(), graph, scratch);
	}
}
#endif
