# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 13:01:30 2020

@author: Kiri Daust
"""

"""Vehicle Routing Problem for PEM sampling"""

from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp
import numpy as np

# file = open("dMatSave.csv")
# dat = np.loadtxt(file, delimiter=",")
# num_days = int(10)
# start = [25,25,25,25,25,25,25,25,25,25]
# end = start
# max_cost = int(480)
# plot_time = int(2)
# GSC = int(10)
# arbDepot = False
# penalty = [29276, 23898, 38238, 30172, 38238, 18521, 28379, 19417, 19417, 
# 12248, 19417, 14040, 14040, 14040, 14040, 23898, 9559, 9559, 
# 19417, 20314, 24795, 18521, 29276, 19417, 19417]


def py_mTSP(dat, num_days, start, end, max_cost, plot_time, penalty, arbDepot, GSC = 10):
    dat2 = dat.copy()
    data = {}
    data['distance_matrix'] = dat2
    data['num_vehicles'] = num_days
    data['starts'] = start
    data['ends'] = end
    
    def start_route(num_days):
        # Create the routing index manager.
        manager = pywrapcp.RoutingIndexManager(len(data['distance_matrix']),num_days, data['starts'][0])
        # Create Routing Model.
        routing = pywrapcp.RoutingModel(manager)
        transit_callback_index = routing.RegisterTransitCallback(distance_callback)

        # Define cost of each arc.
        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)
    
        # Add Distance constraint.
        dimension_name = 'Distance'
        routing.AddDimension(
            transit_callback_index,
            0,  # no slack
            max_cost,  # vehicle maximum travel distance
            True,  # start cumul to zero
            dimension_name)
            
        # Allow to drop nodes.
        for node in range(len(penalty)):
            routing.AddDisjunction([manager.NodeToIndex(node)], penalty[node])
    
        distance_dimension = routing.GetDimensionOrDie(dimension_name)
    
        distance_dimension.SetGlobalSpanCostCoefficient(GSC)
    
        search_parameters = pywrapcp.DefaultRoutingSearchParameters()
        search_parameters.time_limit.seconds = 30
        # Setting first solution heuristic.
        search_parameters = pywrapcp.DefaultRoutingSearchParameters()
        search_parameters.first_solution_strategy = (
            routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)
        return(routing, search_parameters, manager)

    # Create and register a transit callback.
    def distance_callback(from_index, to_index):
        """Returns the distance between the two nodes."""
        # Convert from routing variable Index to distance matrix NodeIndex.
        from_node = manager.IndexToNode(from_index)
        to_node = manager.IndexToNode(to_index)
        if(to_node in end): ##need to fix this part
            add = 0
        else:
            add = plot_time
        return data['distance_matrix'][from_node][to_node]+add
      
    def numEmpty(num_veh):
        empty = 0
        for vehicle_id in range(num_veh):
            index = routing.Start(vehicle_id)
            index = solution.Value(routing.NextVar(index))
            if(routing.IsEnd(index)):
                empty = empty+1
        return(empty)
    
    def num_nodes():
        num = 0
        for vehicle_id in range(num_days):
            index = solution.Value(routing.NextVar(routing.Start(vehicle_id)))
            while not routing.IsEnd(index):
                index = solution.Value(routing.NextVar(index))
                num+=1
        return num

    routing, search_parameters, manager = start_route(num_days)
    # Solve the problem.
    solution = routing.SolveWithParameters(search_parameters)
    
    while(numEmpty(num_days) > 0):
      print("reducing days \n")
      num_days = num_days - numEmpty(num_days)
      routing, search_parameters, manager = start_route(num_days)
      solution = routing.SolveWithParameters(search_parameters)
      
    while(num_days > 1 and (num_nodes()/len(penalty)) > 0.95):
        print("final reduction")
        num_days = num_days - 1
        routing, search_parameters, manager = start_route(num_days)
        solution = routing.SolveWithParameters(search_parameters)
    
    ##collect solution and return
    plan_output = dict.fromkeys(range(num_days),None)
    dist_output = dict.fromkeys(range(num_days),None)
    for vehicle_id in range(num_days):
        index = routing.Start(vehicle_id)
        temp = []
        rd = 0
        while not routing.IsEnd(index):
            temp.append(manager.IndexToNode(index))
            previous_index = index
            index = solution.Value(routing.NextVar(index))
            rd += routing.GetArcCostForVehicle(
                previous_index, index, vehicle_id)
        temp.append(manager.IndexToNode(index))
        if(arbDepot):
            temp = temp[:-1]
        plan_output[vehicle_id] = temp
        dist_output[vehicle_id] = rd
    
    return(plan_output,dist_output,num_days)
