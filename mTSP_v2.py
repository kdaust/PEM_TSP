# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 13:01:30 2020

@author: kirid
"""

"""Vehicles Routing Problem with Delay"""

from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp

def py_mTSP2(dat, num_days, start_node, max_cost, plot_time):
    data = {}
    data['distance_matrix'] = dat
    data['num_vehicles'] = num_days
    data['depot'] = start_node
    data['sample_time'] = plot_time
    
    # Create the routing index manager.
    manager = pywrapcp.RoutingIndexManager(len(data['distance_matrix']),data['num_vehicles'], data['depot'])
    # Create Routing Model.
    routing = pywrapcp.RoutingModel(manager)

    # Create and register a transit callback.
    def distance_callback(from_index, to_index):
        """Returns the distance between the two nodes."""
        # Convert from routing variable Index to distance matrix NodeIndex.
        from_node = manager.IndexToNode(from_index)
        to_node = manager.IndexToNode(to_index)
        return data['distance_matrix'][from_node][to_node]

    transit_callback_index = routing.RegisterTransitCallback(distance_callback)

    # Define cost of each arc.
    routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

    # Add Distance constraint.
    dimension_name = 'Distance'
    routing.AddDimension(
        transit_callback_index,
        5,  # no slack
        max_cost,  # vehicle maximum travel distance
        True,  # start cumul to zero
        dimension_name)
    time_dimension = routing.GetDimensionOrDie(dimension_name)
    
    solver = routing.solver()
    intervals = []
    for i in range(data['num_vehicles']):
        # Add time windows at end of routes.
        intervals.append(
            solver.FixedDurationIntervalVar(
                time_dimension.CumulVar(routing.End(i)),
                data['sample_time'], 'depot_interval'))


    # Setting first solution heuristic.
    search_parameters = pywrapcp.DefaultRoutingSearchParameters()
    search_parameters.first_solution_strategy = (
        routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)

    # Solve the problem.
    solution = routing.SolveWithParameters(search_parameters)
    
    ##collect solution and return
    plan_output = dict.fromkeys(range(data['num_vehicles']),None)
    for vehicle_id in range(data['num_vehicles']):
        index = routing.Start(vehicle_id)
        temp = []
        while not routing.IsEnd(index):
            temp.append(manager.IndexToNode(index))
            index = solution.Value(routing.NextVar(index))
        temp.append(manager.IndexToNode(index))
        plan_output[vehicle_id] = temp
    
    return(plan_output)
