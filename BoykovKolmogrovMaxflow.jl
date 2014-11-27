#= This a translation of Jan Funke's Java implementation of "An Experimental Comparison of
Min-Cut/Max-Flow Algorithms for Energy  Minimization in Vision." Yuri Boykov and Vladimir
Kolmogorov. In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), September 2004 =#

using DataStructures
abstract AbstractNode
abstract AbstractEdge

type NullNode{}<:AbstractNode
    id::Int
    null::Bool
    NullNode(id, null) = new(id, null)

    function NullNode(id)
        return new(id,true)
    end

    function NullNode(id,null)
        return new(id,null)
    end
end


type NullEdge{}<:AbstractEdge
    id::Int
    null::Bool
    NullEdge(id, null) = new(id, null)

    function NullEdge(id)
        return new(id,true)
    end

    function NullEdge(id,null)
        return new(id,null)
    end
end

const NULL_NODE = NullNode(-1)
const NULL_EDGE = NullEdge(-2)

const BACKGROUND = 0
const FOREGROUND = 1

type Node{}<:AbstractNode
    firstOutgoing::AbstractEdge
    parent::AbstractEdge
    next::AbstractNode
    timestamp::Int
    distance::Int
    inSink::Bool
    marked::Bool
    inChangeList::Bool
    residualCapacity::Float64
    id::Int
    null::Bool

    Node(firstOutgoing,parent,next,timestamp,distance,inSink,marked,inChangeList,residualCapacity, id, null) = new(firstOutgoing,parent,next,timestamp,distance,inSink,marked,inChangeList,residualCapacity, id, null)

    function Node(id::Int)
        return new(NULL_EDGE,NULL_EDGE,NULL_NODE,0,0,false,false,false,0.0,id,false)
    end
end

type Edge{}<:AbstractEdge
    head::AbstractNode
    next::AbstractEdge
    sister::AbstractEdge
    residualCapacity::Float64
    id::Int
    null::Bool

    Edge(head,next,sister,residualCapacity,id,null) = new(head,next,sister,residualCapacity,id,null)
    function Edge(id::Int)
        return new(NULL_NODE,NULL_EDGE,NULL_EDGE,0.0,id,false)
    end

end

const TERMINAL_EDGE = Edge(-3)
const ORPHAN_EDGE = Edge(-4)

type GraphCut
    numNodes::Int
    numEdges::Int
    nodes::Array{AbstractNode}
    edges::Array{AbstractEdge}
    edgeNum::Int
    totalFlow::Float64
    maxflowIteration::Int
    activeQueueFirst::Array{AbstractNode}
    activeQueueLast::Array{AbstractNode}
    orphans::Array{AbstractNode}
    time::Int

    GraphCut(numNodes,numEdges,nodes,edges,edgeNum,totalFlow,maxflowIteration,activeQueueFirst,activeQueueLast,orphans,time) = new(numNodes,numEdges,nodes,edges,edgeNum,totalFlow,maxflowIteration,activeQueueFirst,activeQueueLast,orphans,time)

    function GraphCut(numNodes::Int,numEdges::Int)
        nodes = Array(AbstractNode, numNodes)
        for i=1:numNodes
            nodes[i] = Node(i)
        end
        edges = Array(AbstractEdge, 2*numEdges)
        edgeNum = 1
        totalFlow = 0.0
        maxflowIteration = 1
        activeQueueFirst = Array(AbstractNode,2)
        for i=1:2
            activeQueueFirst[i] = NULL_NODE
        end
        activeQueueLast = Array(AbstractNode,2)
        for i=1:2
            activeQueueLast[i] = NULL_NODE
        end
        orphans = AbstractNode[]
        time = 0

        return new(numNodes,numEdges,nodes,edges,edgeNum,totalFlow,maxflowIteration,activeQueueFirst,activeQueueLast,orphans,time)
    end
end

function setTerminalWeights(graph::GraphCut, nodeId::Int, source::Float64, sink::Float64)
    assert (nodeId >= 1 && nodeId <= graph.numNodes)
    delta = graph.nodes[nodeId].residualCapacity
    if delta > 0.0
        source += delta
    else
        source -= delta
    end

    graph.totalFlow += (source < sink) ? source : sink
    graph.nodes[nodeId].residualCapacity = source - sink
end

function setEdgeWeight(graph::GraphCut, nodeId1::Int, nodeId2::Int, weight1to2::Float64, weight2to1::Float64)
    assert (nodeId1 >= 1 && nodeId1 <= graph.numNodes)
    assert (nodeId2 >= 1 && nodeId2 <= graph.numNodes)
    assert (nodeId1 != nodeId2)
    assert (weight1to2 >= 0.0)
    assert (weight2to1 >= 0.0)
    #assert (graph.edgeNum <= graph.numEdges - 2)

    edge = Edge(graph.edgeNum)
    graph.edges[graph.edgeNum] = edge
    graph.edgeNum += 1
    reverseEdge = Edge(graph.edgeNum)
    graph.edges[graph.edgeNum] = reverseEdge
    graph.edgeNum += 1

    node1 = graph.nodes[nodeId1]
    node2 = graph.nodes[nodeId2]

    edge.sister = reverseEdge
    reverseEdge.sister = edge

    edge.next = node1.firstOutgoing
    node1.firstOutgoing = edge

    reverseEdge.next = node2.firstOutgoing
    node2.firstOutgoing = reverseEdge

    edge.head = node2
    reverseEdge.head = node1

    edge.residualCapacity = weight1to2
    reverseEdge.residualCapacity = weight2to1
end

function setEdgeWeight(graph::GraphCut, nodeId1::Int, nodeId2::Int, weight::Float64)
    setEdgeWeight(graph, nodeId1, nodeId2, weight, weight)
end

function getTerminal(graph::GraphCut, nodeId::Int)
    assert (nodeId >= 1 && nodeId <= graph.numNodes);

    node = graph.nodes[nodeId]

    if !node.parent.null
       return node.inSink ? BACKGROUND : FOREGROUND
    else
       return BACKGROUND
    end
end

function markNode(graph::GraphCut, nodeId::Int)
    assert (nodeId >= 1 && nodeId <= numNodes)
    node = graph.nodes[nodeId]

    if node.next.null
        if !graph.activeQueueLast[2].null
            graph.activeQueueLast[2].next = node
        else
            graph.activeQueueFirst[2] = node
        end

        graph.activeQueueLast[2] = node;
        node.next = node
    end
    node.marked = true
end

function setNodeActive(graph::GraphCut, node::Node)
    if node.next.null
        if !graph.activeQueueLast[2].null
            graph.activeQueueLast[2].next = node
        else
            graph.activeQueueFirst[2] = node
        end

        graph.activeQueueLast[2] = node;
        node.next = node
    end
end

function getNextActiveNode(graph::GraphCut)
    while true
        node = graph.activeQueueFirst[1]

        if node.null
            node = graph.activeQueueFirst[2]

            graph.activeQueueFirst[1] = graph.activeQueueFirst[2]
            graph.activeQueueLast[1] = graph.activeQueueLast[2]
            graph.activeQueueFirst[2] = NULL_NODE
            graph.activeQueueLast[2] = NULL_NODE

            if node.null
                return NULL_NODE
            end
        end

        if node.next.id == node.id
            graph.activeQueueFirst[1] = NULL_NODE
            graph.activeQueueLast[1] = NULL_NODE
        else
            graph.activeQueueFirst[1] = node.next
        end

        node.next = NULL_NODE

        if !node.parent.null
            return node;
        end
    end
end

function addOrphanAtFront(graph::GraphCut, node::Node)
    node.parent = ORPHAN_EDGE
    unshift!(graph.orphans, node) # add to front
end

function addOrphanAtBack(graph::GraphCut, node::Node)
    node.parent = ORPHAN_EDGE
    push!(graph.orphans, node) # add to back
end

function addToChangedList(graph::GraphCut, node::Node)
    node.inChangeList = true
end

function maxflowInit(graph::GraphCut)
    graph.activeQueueFirst[1] = NULL_NODE
    graph.activeQueueLast[1] = NULL_NODE
    graph.activeQueueFirst[2] = NULL_NODE
    graph.activeQueueLast[2] = NULL_NODE

    graph.orphans = AbstractNode[]
    graph.time = 0

    for node in graph.nodes
        node.next = NULL_NODE
        node.marked = false
        node.inChangeList = false
        node.timestamp = graph.time

        if node.residualCapacity > 0.0
            node.inSink = false
            node.parent = TERMINAL_EDGE
            setNodeActive(graph, node)
            node.distance = 1
        elseif node.residualCapacity < 0.0
            node.inSink = true
            node.parent = TERMINAL_EDGE
            setNodeActive(graph, node)
            node.distance = 1
        else
            node.parent = NULL_EDGE
        end
    end
end

function maxflowReuseTreesInit(graph::GraphCut)
        Node::node1
        Node::node2

        queueStart = graph.activeQueueFirst[2]

        Edge::edge

        graph.activeQueueFirst[1] = NULL_NODE
        graph.activeQueueLast[1] = NULL_NODE
        graph.activeQueueFirst[2] = NULL_NODE
        graph.activeQueueLast[2] = NULL_NODE

        graph.orphans = AbstractNode[]

        graph.time += 1

        while !(node1 = queueStart).null
            queueStart = node1.next

            if queueStart.id == node1.id
                queueStart = NULL_NODE
            end

            node1.next = NULL_NODE
            node1.marked = false
            setNodeActive(graph, node1);

            if node1.getResidualCapacity() == 0.0
                if !node1.parent.null
                    addOrphanAtBack(graph, node1)
                end
                continue
            end

            if node1.residualCapacity > 0.0
                if node1.parent.null || node1.inSink
                    node1.inSink = false
                    edge = node1.firstOutgoing
                    while !edge.null
                        node2 = edge.head
                        if node2.marked
                            if node2.parent.id == edge.sister.id
                                addOrphanAtBack(graph, node2)
                            end

                            if !node2.parent.null && node2.inSink && edge.residualCapacity > 0.0
                                setNodeActive(graph, node2)
                            end

                        end
                        edge = edge.next
                    end
                    addToChangedList(graph, node1)
            else
                if node1.parent.null || !node.inSink
                   node1.inSink = true
                   edge = node1.firstOutgoing
                   while !edge.null
                        node2 = edge.head
                        if !node2.marked
                            if node2.parent.id == edge.sister.id
                                addOrphanAtBack(graph, node2)
                            end
                            if !node2.parent.null && !node2.inSink && edge.sister.residualCapacity > 0
                               setActiveNode(graph, node2)
                            end
                        end
                        edge = edge.next
                   end
                   addToChangedList(graph, node1)
                end
            end
            node1.parent = TERMINAL_EDGE
            node1.timestamp = time
            node1.distance = 1.0
        end


        while length(graph.orphans) > 0
            orphan = shift!(graph.orphans)
            if orphan.inSink
                processSinkOrphan(graph, orphan)
            else
                processSourceOrphan(graph, orphan)
            end
        end
    end
end


function augment(graph::GraphCut, middle::Edge)
    bottleneck = middle.residualCapacity
    node = middle.sister.head
    while true
        edge = node.parent

        if edge.id==TERMINAL_EDGE.id
            break
        end

        if bottleneck > edge.sister.residualCapacity
            bottleneck = edge.sister.residualCapacity
        end

        node = edge.head
    end

    if bottleneck > node.residualCapacity
       bottleneck = node.residualCapacity
    end

    node = middle.head
    while true
        edge = node.parent

        if edge.id == TERMINAL_EDGE.id
            break
        end

        if bottleneck > edge.residualCapacity
            bottleneck = edge.residualCapacity
        end

        node = edge.head
    end

    if bottleneck > -node.residualCapacity
        bottleneck = -node.residualCapacity
    end

    middle.sister.residualCapacity += bottleneck
    middle.residualCapacity -= bottleneck

    node = middle.sister.head
    while true
       edge = node.parent

       if edge.id == TERMINAL_EDGE.id
           break
       end

       edge.residualCapacity += bottleneck
       edge.sister.residualCapacity -= bottleneck
       if edge.sister.residualCapacity == 0.0
            addOrphanAtFront(graph, node)
       end

       node = edge.head
    end
    node.residualCapacity -= bottleneck
    if node.residualCapacity == 0.0
        addOrphanAtFront(graph, node)
    end

    node = middle.head
    while true
       edge = node.parent
       if edge.id == TERMINAL_EDGE.id
            break
       end
       edge.sister.residualCapacity += bottleneck
       edge.residualCapacity -= bottleneck
       if edge.residualCapacity == 0.0
            addOrphanAtFront(graph, node)
       end
       node = edge.head
    end

    node.residualCapacity += bottleneck

    if node.residualCapacity == 0.0
        addOrphanAtFront(graph, node)
    end

    graph.totalFlow += bottleneck
end

function processSourceOrphan(graph::GraphCut, orphan::Node)

    bestEdge = NULL_EDGE
    minDistance = typemax(Int32)

    orphanEdge = orphan.firstOutgoing
    while !orphanEdge.null
        if orphanEdge.sister.residualCapacity != 0.0
            node = orphanEdge.head
            parentEdge = node.parent
            if !node.inSink && !parentEdge.null
                distance = 0
                while true
                    if node.timestamp == graph.time
                        distance += node.distance
                        break
                    end
                    parentEdge = node.parent
                    distance += 1
                    if parentEdge.id == TERMINAL_EDGE.id
                       node.timestamp = graph.time
                       node.distance = 1
                       break;
                    end

                    if parentEdge.id == ORPHAN_EDGE.id
                       distance = typemax(Int32)
                       break;
                    end
                    node = parentEdge.head
                end

                if distance < typemax(Int32)
                    if distance < minDistance
                        bestEdge = orphanEdge
                        minDistance = distance
                    end

                    node = orphanEdge.head
                    while node.timestamp != graph.time
                        node.timestamp = graph.time
                        node.distance = distance
                        distance -= 1
                        node = node.parent.head
                    end
                end
            end
        end
        orphanEdge = orphanEdge.next
    end

    orphan.parent = bestEdge
    if !bestEdge.null
       orphan.timestamp = graph.time
       orphan.distance = minDistance + 1
    else
        addToChangedList(graph, orphan)

        orphanEdge = orphan.firstOutgoing
        while !orphanEdge.null
            node = orphanEdge.head
            parentEdge = node.parent

            if !node.inSink && !parentEdge.null
                if orphanEdge.sister.residualCapacity != 0.0
                    setNodeActive(graph, node)
                end

                if parentEdge.id != TERMINAL_EDGE.id && parentEdge.id != ORPHAN_EDGE.id && parentEdge.head.id == orphan.id
                    addOrphanAtBack(graph, node)
                end
            end

            orphanEdge = orphanEdge.next
        end
    end
end

function processSinkOrphan(graph, orphan::Node)
    bestEdge = NULL_EDGE
    minDistance = typemax(Int32)

    orphanEdge = orphan.firstOutgoing
    while !orphanEdge.null
        if orphanEdge.residualCapacity != 0.0
           node = orphanEdge.head
           parentEdge = node.parent

           if node.inSink && !parentEdge.null
                distance = 0
                while true
                    if node.timestamp == graph.time
                        distance += node.distance
                        break
                    end
                    parentEdge = node.parent
                    distance += 1
                    if parentEdge.id == TERMINAL_EDGE.id
                        node.timestamp = graph.time
                        node.distance = 1
                        break
                    end
                    if parentEdge.id == ORPHAN_EDGE.id
                       distance = typemax(Int32)
                       break
                    end
                    node = parentEdge.head
                end

                if distance < typemax(Int32)
                    if distance < minDistance
                        bestEdge = orphanEdge
                        minDistance = distance
                    end

                    node = orphanEdge.head
                    while node.timestamp != graph.time
                        node.timestamp = graph.time
                        node.distance = distance
                        distance -= 1
                        node = node.parent.head
                    end
                end
           end
        end
        orphanEdge = orphanEdge.next
    end

    orphan.parent = bestEdge
    if !bestEdge.null
        orphan.timestamp = graph.time
        orphan.distance = minDistance + 1
    else
        addToChangedList(graph, orphan)

        orphanEdge = orphan.firstOutgoing
        while !orphanEdge.null
            node = orphanEdge.head
            parentEdge = node.parent
            if node.inSink && !parentEdge.null
                if orphanEdge.residualCapacity != 0
                    setNodeActive(graph, node)
                end
                if parentEdge.id != TERMINAL_EDGE.id && parentEdge.id != ORPHAN_EDGE.id && parentEdge.head.id == orphan.id
                    addOrphanAtBack(graph, node)
                end
            end

            orphanEdge = orphanEdge.next
        end
    end
end

function computeMaximumFlow(graph::GraphCut, reuseTrees::Bool, changedNodes::Array{Int})
    if graph.maxflowIteration == 1
        reuseTrees = false
    end

    if reuseTrees
        maxflowReuseTreesInit(graph)
    else
        maxflowInit(graph)
    end

    currentNode = NULL_NODE
    edge = NULL_EDGE
    while true
        activeNode = currentNode

        if !activeNode.null
            activeNode.next = NULL_NODE
            if activeNode.parent.null
                activeNode = NULL_NODE
            end
        end
        if activeNode.null
            activeNode = getNextActiveNode(graph)
            if activeNode.null
                break
            end
        end

        if !activeNode.inSink
            edge = activeNode.firstOutgoing
            while !edge.null
                if edge.residualCapacity != 0.0
                    headNode = edge.head

                    if headNode.parent.null
                        headNode.inSink = false
                        headNode.parent = edge.sister
                        headNode.timestamp = activeNode.timestamp
                        headNode.distance = activeNode.distance + 1
                        setNodeActive(graph, headNode)
                        addToChangedList(graph, headNode)
                    elseif headNode.inSink
                        break;
                    elseif headNode.timestamp <= activeNode.timestamp && headNode.distance > activeNode.distance
                        headNode.parent = edge.sister
                        headNode.timestamp = activeNode.timestamp
                        headNode.distance = activeNode.distance + 1
                    end
                end
                edge = edge.next
            end
        else
            edge = activeNode.firstOutgoing
            while !edge.null
                if edge.sister.residualCapacity != 0.0
                    headNode = edge.head

                    if headNode.parent.null
                        headNode.inSink = true
                        headNode.parent = edge.sister
                        headNode.timestamp = activeNode.timestamp
                        headNode.distance = activeNode.distance + 1
                        setNodeActive(graph, headNode)
                        addToChangedList(graph, headNode)
                    elseif !headNode.inSink
                        edge = edge.sister
                        break;
                    elseif headNode.timestamp <= activeNode.timestamp && headNode.distance > activeNode.distance
                        headNode.parent = edge.sister
                        headNode.timestamp = activeNode.timestamp
                        headNode.distance = activeNode.distance + 1
                    end
                end
                edge = edge.next
            end
        end

        graph.time += 1

        if !edge.null
            activeNode.next = activeNode
            currentNode = activeNode
            augment(graph, edge)

            while length(graph.orphans) > 0
                orphan = shift!(graph.orphans) #remove front
                if orphan.inSink
                     processSinkOrphan(graph, orphan)
                else
                     processSourceOrphan(graph, orphan)
                end
            end
        else
            currentNode = NULL_NODE
        end
    end

    graph.maxflowIteration += 1

    if length(changedNodes)>=0
         changedNodes = Int[]
         for i=1:length(graph.nodes)
            if graph.nodes[i].inChangeList
                push!(changedNodes,i)
            end
         end
    end

    return graph.totalFlow
end

