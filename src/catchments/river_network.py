import os
import geopandas as gpd
import networkx as nx
from shapely.geometry import Point, LineString
from shapely.ops import split
from collections import defaultdict


class RiverNetwork:
    """
    A class to manage river network analysis with stations.
    
    Attributes:
    -----------
    streams_gdf : gpd.GeoDataFrame
        Original stream network
    stations_gdf : gpd.GeoDataFrame
        Station points
    streams_split : gpd.GeoDataFrame
        Stream network split at station locations
    graph : nx.DiGraph
        Directed graph representation of the river network
    station_locations : dict
        Mapping of station points to station IDs
    """
    
    def __init__(self, streams_gdf: gpd.GeoDataFrame, stations_gdf: gpd.GeoDataFrame, 
                 strahler_col: str = 'strord'):
        """
        Initialize RiverNetwork with streams and stations.
        
        Parameters:
        -----------
        streams_gdf : gpd.GeoDataFrame
            River network with LineString geometries
        stations_gdf : gpd.GeoDataFrame
            Station points to split the network
        strahler_col : str
            Column name containing Strahler order (default: 'strord')
        """
        # Ensure same CRS
        if streams_gdf.crs != stations_gdf.crs:
            raise ValueError("Streams and stations must have the same CRS")
        
        self.streams_gdf = streams_gdf.copy()
        self.stations_gdf = stations_gdf.copy()
        self.strahler_col = strahler_col
        self.streams_split = None
        self.graph = None
        
        # Store station locations for later node attribution
        self.station_locations = {}
        for _, station in stations_gdf.iterrows():
            point = station.geometry
            # Round coordinates to match node IDs
            coord_key = (round(point.x, 6), round(point.y, 6))
            self.station_locations[coord_key] = station['codigoestacao']
        
        print(f"RiverNetwork initialized with {len(streams_gdf)} streams and {len(stations_gdf)} stations")
    
    
    def split_streams_by_stations(self, tolerance: float = 10, 
                            save_data: bool = False,
                            output_dir: str = None,
                            file_name: str = 'streams_split_by_station'):
        """
        Split stream geometries wherever a station intersects them.
        
        Parameters:
        -----------
        tolerance : float
            Distance tolerance for snapping stations to streams (in CRS units)
        save_data : bool
            Whether to save the result to file
        output_dir : str
            Directory to save output (default: './outputs')
        file_name : str
            Name of output file
            
        Returns:
        --------
        gpd.GeoDataFrame
            Stream network split at station locations
        """
        # Setup output directory
        if save_data:
            if output_dir is None:
                output_dir = os.path.join(os.getcwd(), 'outputs')
            os.makedirs(output_dir, exist_ok=True)
        
        out_rows = []
        split_count = 0
        no_split_count = 0
        
        # Build spatial index for stations
        print("Building spatial index for stations.")
        stations_sindex = self.stations_gdf.sindex
        
        print("Splitting network on point intersection.")
        for _, stream in self.streams_gdf.iterrows():
            geom = stream.geometry
            
            # Find candidate stations
            possible_idx = list(stations_sindex.intersection(geom.bounds))
            
            if not possible_idx:
                # No stations nearby, keep original
                row = stream.to_dict()
                out_rows.append(row)
                continue
            
            # Get stations that are within tolerance and snap them ONTO the line
            stations = self.stations_gdf.iloc[possible_idx]
            intersecting_stations = []
            
            for _, s in stations.iterrows():
                station_point = s.geometry
                # Check distance to line
                dist = geom.distance(station_point)
                
                if dist <= tolerance:
                    # Snap the point ONTO the line using projection
                    # project() gives distance along line, interpolate() gets point at that distance
                    snapped_point = geom.interpolate(geom.project(station_point))
                    intersecting_stations.append(snapped_point)
            
            if not intersecting_stations:
                # No stations within tolerance, keep original
                row = stream.to_dict()
                out_rows.append(row)
                no_split_count += 1
                continue
            
            # Split by all stations at once
            try:
                from shapely.geometry import MultiPoint
                station_points = MultiPoint(intersecting_stations)
                parts = split(geom, station_points)
                
                # Check if split actually worked
                if len(parts.geoms) == 1:
                    # No split occurred (shouldn't happen but just in case)
                    print(f"Warning: Split failed for stream, keeping original")
                    row = stream.to_dict()
                    out_rows.append(row)
                    no_split_count += 1
                else:
                    # Split successful - add all parts
                    for part in parts.geoms:
                        row = stream.drop("geometry").to_dict()
                        row["geometry"] = part
                        out_rows.append(row)
                    split_count += 1
                    
            except Exception as e:
                print(f"Failed to split stream: {e}")
                # Keep original geometry
                row = stream.to_dict()
                out_rows.append(row)
                no_split_count += 1
        
        self.streams_split = gpd.GeoDataFrame(out_rows, crs=self.streams_gdf.crs)
        
        # Write dataset to file
        if save_data:
            try:
                print(f"Saving split streams to {output_dir}")
                self.streams_split.to_file(
                    os.path.join(output_dir, f'{file_name}.geojson'), 
                    driver='GeoJSON'
                )
            except Exception as e:
                print(f'File could not be saved: {e}')
        
        print(f"  Streams that were split: {split_count}")
        
        return self.streams_split
    
    
    def fix_river_direction(self, gdf=None):
        """
        Automatically orient river segments so they flow from low to high Strahler order.
        Rivers flow from tributaries (low order) to main stem (high order).
        
        Parameters:
        -----------
        gdf : gpd.GeoDataFrame, optional
            GeoDataFrame to fix. If None, uses streams_split or streams_gdf
            
        Returns:
        --------
        gpd.GeoDataFrame
            GeoDataFrame with corrected flow directions
        """
        if gdf is None:
            gdf = self.streams_split if self.streams_split is not None else self.streams_gdf
        
        gdf = gdf.copy()
        
        # Build connectivity map
        endpoint_to_segments = defaultdict(list)
        
        for idx, row in gdf.iterrows():
            try:
                geom = row.geometry
                start = (round(geom.coords[0][0], 6), round(geom.coords[0][1], 6))
                end = (round(geom.coords[-1][0], 6), round(geom.coords[-1][1], 6))
                
                endpoint_to_segments[start].append((idx, 'start'))
                endpoint_to_segments[end].append((idx, 'end'))
            except Exception as e:
                print(f'Failed to get coordinates for geometry {row.geometry}, {e}')
        
        # Fix orientation based on Strahler order
        for point, connections in endpoint_to_segments.items():
            if len(connections) < 2:
                continue
            
            # Get segments and their Strahler orders
            segments_info = []
            for seg_idx, endpoint_type in connections:
                strahler = gdf.loc[seg_idx, self.strahler_col]
                segments_info.append((seg_idx, endpoint_type, strahler))
            
            # Find if there's a clear high-order segment
            strahler_orders = [s[2] for s in segments_info]
            max_strahler = max(strahler_orders)
            
            # If there's a dominant high-order stream, orient others toward it
            low_order_segments = [s for s in segments_info if s[2] < max_strahler]
            
            # Low order segments should flow INTO (end at) the junction
            for seg_idx, endpoint_type, _ in low_order_segments:
                if endpoint_type == 'start':
                    # Reverse this segment so it ends at the junction
                    gdf.at[seg_idx, 'geometry'] = LineString(
                        list(gdf.loc[seg_idx, 'geometry'].coords)[::-1]
                    )
        
        # Update the appropriate attribute
        if self.streams_split is not None:
            self.streams_split = gdf
        else:
            self.streams_gdf = gdf
        
        print("River directions corrected based on Strahler order")
        return gdf
    
    
    def create_graph(self, direction='downstream', fix_direction=True):
        """
        Create a directed graph from the river network.
        Nodes at station locations will be marked with station_id.
        
        Parameters:
        -----------
        direction : str
            'downstream' or 'upstream' for edge direction
        fix_direction : bool
            Whether to auto-correct flow directions before creating graph
            
        Returns:
        --------
        nx.DiGraph
            Directed graph with river segments
        """
        # Use split streams if available, otherwise use original
        gdf = self.streams_split if self.streams_split is not None else self.streams_gdf
        
        # Fix directions if requested
        if fix_direction:
            gdf = self.fix_river_direction(gdf)
        
        G = nx.DiGraph()
        
        # First pass: Create all nodes and edges
        for idx, row in gdf.iterrows():
            try:
                geom = row.geometry
                strahler = row[self.strahler_col]
                
                # Get start and end points
                start_point = Point(geom.coords[0])
                end_point = Point(geom.coords[-1])
                
                # Create node IDs from coordinates
                start_id = (round(start_point.x, 6), round(start_point.y, 6))
                end_id = (round(end_point.x, 6), round(end_point.y, 6))
                
                # Add nodes with basic attributes
                if not G.has_node(start_id):
                    G.add_node(start_id, pos=start_point, type='junction')
                if not G.has_node(end_id):
                    G.add_node(end_id, pos=end_point, type='junction')
                
                # Add edge with attributes (no station_id here)
                edge_attrs = {
                    'segment_id': idx,
                    'strahler': strahler,
                    'geometry': geom,
                    'length': geom.length
                }
                
                if direction == 'downstream':
                    G.add_edge(start_id, end_id, **edge_attrs)
                else:  # upstream
                    G.add_edge(end_id, start_id, **edge_attrs)
            
            except Exception as e:
                print(f'Failed to process geometry at index {idx}: {e}')
                continue
        
        # Second pass: Mark nodes that are stations
        stations_found = 0
        for node_id in G.nodes():
            if node_id in self.station_locations:
                G.nodes[node_id]['type'] = 'station'
                G.nodes[node_id]['station_id'] = self.station_locations[node_id]
                stations_found += 1
        
        self.graph = G
        print(f"Graph created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
        print(f"Stations found in graph: {stations_found}/{len(self.station_locations)}")
        return G
    
    
    def get_network_stats(self):
        """Get summary statistics about the river network."""
        stats = {
            'n_streams': len(self.streams_gdf),
            'n_stations': len(self.stations_gdf),
        }
        
        if self.streams_split is not None:
            stats['n_streams_splitted'] = len(self.streams_split)
        
        if self.graph is not None:
            stats['n_nodes'] = self.graph.number_of_nodes()
            stats['n_edges'] = self.graph.number_of_edges()
            
            # Find outlets
            outlets = [n for n in self.graph.nodes() if self.graph.out_degree(n) == 0]
            stats['n_outlets'] = len(outlets)
            
            # Find headwaters
            headwaters = [n for n in self.graph.nodes() if self.graph.in_degree(n) == 0]
            stats['n_headwaters'] = len(headwaters)
            
            # Count station nodes (unique station_id values)
            station_nodes = [n for n, d in self.graph.nodes(data=True) 
                           if d.get('type') == 'station']
            stats['n_station_nodes'] = len(station_nodes)
        
        return stats
    

    def get_station_nodes(self):
        """
        Get all nodes that represent stations.
        
        Returns:
        --------
        dict
            Dictionary mapping node coordinates to station_id
        """
        if self.graph is None:
            raise ValueError("Graph not created. Call create_graph() first.")
        
        station_nodes = {}
        for node, data in self.graph.nodes(data=True):
            if data.get('type') == 'station' and 'station_id' in data:
                station_nodes[node] = data['station_id']
        
        return station_nodes
    
    
    def list_all_station_ids(self):
        """
        List all unique station IDs in the graph.
        
        Returns:
        --------
        list
            List of unique station IDs
        """
        if self.graph is None:
            raise ValueError("Graph not created. Call create_graph() first.")
        
        station_ids = []
        for node, data in self.graph.nodes(data=True):
            if 'station_id' in data and data['station_id'] is not None:
                station_ids.append(data['station_id'])
        
        return station_ids  # Already unique since each node can only have one station_id
    
    
    def graph_to_geodataframe(self, include_cumulative_length=True, 
                            start_nodes=None, direction='downstream'):
        """
        Convert the river network graph to a GeoDataFrame with cumulative lengths.
        
        Parameters:
        -----------
        include_cumulative_length : bool
            Whether to calculate cumulative river length from headwaters/outlets
        start_nodes : list, optional
            List of nodes to start cumulative length calculation from.
            If None, uses headwaters (for downstream) or outlets (for upstream)
        direction : str
            'downstream' (cumulative from headwaters) or 'upstream' (cumulative from outlets)
            
        Returns:
        --------
        gpd.GeoDataFrame
            GeoDataFrame with all edges and their attributes including cumulative length
        """
        if self.graph is None:
            raise ValueError("Graph not created. Call create_graph() first.")
        
        edges_data = []
        
        # First, collect all edge data
        for u, v, data in self.graph.edges(data=True):
            edge_dict = {
                'geometry': data['geometry'],
                'from_node': str(u),
                'to_node': str(v),
                'segment_id': data['segment_id'],
                'strahler': data['strahler'],
                'length': data['length'],
                'length_km': data['length'] / 1000
            }
            
            # Add any additional edge attributes
            for key in ['closest_station', 'cumulative_distance']:
                if key in data:
                    edge_dict[key] = data[key]
            
            edges_data.append(edge_dict)
        
        edges_gdf = gpd.GeoDataFrame(edges_data, crs=self.streams_gdf.crs)
        
        # Calculate cumulative length if requested
        if include_cumulative_length:
            edges_gdf = self._add_cumulative_length(edges_gdf, start_nodes, direction)
        
        print(f"Graph converted to GeoDataFrame: {len(edges_gdf)} edges")
        
        return edges_gdf


    def _add_cumulative_length(self, edges_gdf, start_nodes=None, direction='downstream'):
        """
        Helper method to add cumulative length to edges GeoDataFrame.
        
        Parameters:
        -----------
        edges_gdf : gpd.GeoDataFrame
            GeoDataFrame with edges
        start_nodes : list, optional
            Starting nodes for cumulative calculation
        direction : str
            'downstream' or 'upstream'
            
        Returns:
        --------
        gpd.GeoDataFrame
            GeoDataFrame with cumulative_length column added
        """
        # Initialize cumulative length
        edges_gdf['cumulative_length'] = 0.0
        
        # Determine start nodes
        if start_nodes is None:
            if direction == 'downstream':
                # Start from headwaters (nodes with no incoming edges)
                start_nodes = [n for n in self.graph.nodes() 
                            if self.graph.in_degree(n) == 0]
            else:  # upstream
                # Start from outlets (nodes with no outgoing edges)
                start_nodes = [n for n in self.graph.nodes() 
                            if self.graph.out_degree(n) == 0]
        
        print(f"Calculating cumulative length from {len(start_nodes)} start nodes...")
        
        # Track cumulative distance to each node
        node_cumulative_dist = {node: 0.0 for node in start_nodes}
        
        # BFS traversal to calculate cumulative lengths
        visited_edges = set()
        queue = list(start_nodes)
        visited_nodes = set(start_nodes)
        
        while queue:
            current_node = queue.pop(0)
            current_cumulative = node_cumulative_dist[current_node]
            
            # Get next nodes based on direction
            if direction == 'downstream':
                next_edges = [
                    (current_node, succ, self.graph.edges[current_node, succ])
                    for succ in self.graph.successors(current_node)
                ]
            else:  # upstream
                next_edges = [
                    (pred, current_node, self.graph.edges[pred, current_node])
                    for pred in self.graph.predecessors(current_node)
                ]
            
            for u, v, edge_data in next_edges:
                edge_key = (u, v)
                
                if edge_key in visited_edges:
                    continue
                
                # Calculate cumulative length for this edge
                edge_length = edge_data['length'] / 1000  # Convert to km
                new_cumulative = current_cumulative + edge_length
                
                # Update the GeoDataFrame
                mask = (edges_gdf['from_node'] == str(u)) & (edges_gdf['to_node'] == str(v))
                edges_gdf.loc[mask, 'cumulative_length'] = new_cumulative
                
                visited_edges.add(edge_key)
                
                # Update next node's cumulative distance
                next_node = v if direction == 'downstream' else u
                
                if next_node not in node_cumulative_dist:
                    node_cumulative_dist[next_node] = new_cumulative
                else:
                    # Take the maximum cumulative distance if multiple paths lead here
                    node_cumulative_dist[next_node] = max(
                        node_cumulative_dist[next_node], 
                        new_cumulative
                    )
                
                # Add to queue if not visited
                if next_node not in visited_nodes:
                    visited_nodes.add(next_node)
                    queue.append(next_node)
        
        print(f"Cumulative length calculated for {len(visited_edges)} edges")
        
        return edges_gdf


    def get_upstream_subgraph(self, station_id, include_station_node=True):
        """
        Get the complete upstream subgraph for a given station.
        Uses iterative BFS to avoid recursion limits.
        
        Parameters:
        -----------
        station_id : str
            ID of the reference station
        include_station_node : bool
            Whether to include the station node itself in the subgraph
            
        Returns:
        --------
        nx.DiGraph
            Subgraph containing all upstream nodes and edges
        """
        if self.graph is None:
            raise ValueError("Graph not created. Call create_graph() first.")
        
        # Find the node for this station
        station_node = None
        for node, data in self.graph.nodes(data=True):
            if data.get('station_id') == station_id:
                station_node = node
                break
        
        if station_node is None:
            raise ValueError(f"Station {station_id} not found in graph")
        
        # Use iterative BFS to find all upstream nodes (avoiding recursion)
        upstream_nodes = set()
        queue = [station_node]
        visited = {station_node}
        
        while queue:
            current_node = queue.pop(0)
            
            # Get all predecessors (nodes that flow into current_node)
            predecessors = list(self.graph.predecessors(current_node))
            
            for pred in predecessors:
                if pred not in visited:
                    visited.add(pred)
                    upstream_nodes.add(pred)
                    queue.append(pred)
        
        # Include the station node if requested
        if include_station_node:
            subgraph_nodes = upstream_nodes | {station_node}
        else:
            subgraph_nodes = upstream_nodes
        
        # Create subgraph with all upstream nodes and their connecting edges
        upstream_subgraph = self.graph.subgraph(subgraph_nodes).copy()
        
        print(f"Upstream subgraph for station {station_id}:")
        print(f"  Nodes: {upstream_subgraph.number_of_nodes()}")
        print(f"  Edges: {upstream_subgraph.number_of_edges()}")
        
        # Count stations in subgraph
        station_count = sum(1 for n, d in upstream_subgraph.nodes(data=True) 
                        if d.get('type') == 'station')
        print(f"  Stations: {station_count}")
        
        return upstream_subgraph


    def visualize_upstream_subgraph(self, station_id, include_station_node=True, 
                                figsize=(12, 10), save_path=None):
        """
        Visualize the upstream subgraph for a station.
        
        Parameters:
        -----------
        station_id : str
            ID of the reference station
        include_station_node : bool
            Whether to include the station node
        figsize : tuple
            Figure size
        save_path : str, optional
            Path to save figure
        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        import matplotlib.cm as cm
        
        subgraph = self.get_upstream_subgraph(station_id, include_station_node)
        
        # Get positions
        pos = nx.get_node_attributes(subgraph, 'pos')
        pos_coords = {k: (v.x, v.y) for k, v in pos.items()}
        
        # Prepare node colors
        node_colors = []
        node_sizes = []
        for node, data in subgraph.nodes(data=True):
            if data.get('type') == 'station':
                if data.get('station_id') == station_id:
                    node_colors.append('red')  # Target station
                    node_sizes.append(200)
                else:
                    node_colors.append('orange')  # Other stations
                    node_sizes.append(150)
            else:
                node_colors.append('lightblue')  # Regular junctions
                node_sizes.append(50)
        
        # Prepare edge colors based on Strahler order
        strahler_values = [data['strahler'] for u, v, data in subgraph.edges(data=True)]
        norm = Normalize(vmin=min(strahler_values), vmax=max(strahler_values))
        cmap = cm.cool
        edge_colors = [cmap(norm(s)) for s in strahler_values]
        
        # Plot
        fig, ax = plt.subplots(figsize=figsize)
        
        # Draw edges
        nx.draw_networkx_edges(subgraph, pos_coords, edge_color=edge_colors,
                            width=2, alpha=0.6, arrows=True, arrowsize=10,
                            ax=ax)
        
        # Draw nodes
        nx.draw_networkx_nodes(subgraph, pos_coords, node_color=node_colors,
                            node_size=node_sizes, alpha=0.8, ax=ax)
        
        # Add colorbar for Strahler order
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, label='Strahler Order')
        
        ax.set_title(f'Upstream Network for Station {station_id}', fontsize=14, fontweight='bold')
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.axis('equal')
        
        # Add legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='red', 
                markersize=10, label=f'Target Station ({station_id})'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', 
                markersize=10, label='Other Stations'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='lightblue', 
                markersize=7, label='Junctions')
        ]
        ax.legend(handles=legend_elements, loc='best')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {save_path}")
        
        plt.show()



    def get_upstream_geodataframe(self, station_id, include_station_node=True, 
                                include_nodes=False):
        """
        Get upstream subgraph as a GeoDataFrame.
        
        Parameters:
        -----------
        station_id : str
            ID of the reference station
        include_station_node : bool
            Whether to include the station node itself in the subgraph
        include_nodes : bool
            If True, returns (edges_gdf, nodes_gdf). If False, returns only edges_gdf
            
        Returns:
        --------
        gpd.GeoDataFrame or tuple of gpd.GeoDataFrame
            GeoDataFrame(s) with upstream network edges (and optionally nodes)
        """
        # Get the upstream subgraph
        subgraph = self.get_upstream_subgraph(station_id, include_station_node)
        
        # Convert edges to GeoDataFrame
        edges_data = []
        for u, v, data in subgraph.edges(data=True):
            edge_dict = {
                'geometry': data['geometry'],
                'from_node': str(u),
                'to_node': str(v),
                'segment_id': data['segment_id'],
                'strahler': data['strahler'],
                'length': data['length']
            }
            edges_data.append(edge_dict)
        
        edges_gdf = gpd.GeoDataFrame(edges_data, crs=self.streams_gdf.crs)
        
        if not include_nodes:
            print(f"Upstream network for station {station_id}: {len(edges_gdf)} edges")
            return edges_gdf
        
        # Convert nodes to GeoDataFrame
        nodes_data = []
        for node, data in subgraph.nodes(data=True):
            node_dict = {
                'geometry': data['pos'],
                'node_id': str(node),
                'node_type': data.get('type', 'junction'),
                'station_id': data.get('station_id', None),
                'in_degree': subgraph.in_degree(node),
                'out_degree': subgraph.out_degree(node)
            }
            nodes_data.append(node_dict)
        
        nodes_gdf = gpd.GeoDataFrame(nodes_data, crs=self.streams_gdf.crs)
        
        print(f"Upstream network for station {station_id}:")
        print(f"  Edges: {len(edges_gdf)}")
        print(f"  Nodes: {len(nodes_gdf)}")
        
        return edges_gdf, nodes_gdf
    

    def compute_station_influence_zones(self, distance_col='dist_to_station'):
        """
        Compute upstream influence zones for each station by crawling the network.
        
        This method:
        1. Creates a copy of the main graph
        2. For each station, crawls upstream summing lengths
        3. Updates nodes with their closest station and distance
        4. Prioritizes lower Strahler orders at junctions
        5. Stops when cumulative distance exceeds max_distance_km
        
        Parameters:
        -----------
        distance_col : str
            Column name for distance attribute in nodes
            
        Returns:
        --------
        nx.DiGraph
            Annotated graph with station influence zones
        """

        while True:
            try:
                max_distance_km = input(
                    f"Choose the maximum recommended river length distance"
                    f"that should be monitored by a station (in km). Type q to cancel and quit: "
                )

                # Check empty input
                if max_distance_km.lower() == 'q':
                    return None, None
                
                # Check empty input
                if not max_distance_km.strip():
                    raise ValueError("Input cannot be empty.")

                max_distance_km = int(max_distance_km)

                # Check positivity
                if max_distance_km <= 0:
                    raise ValueError("Distance must be a positive integer.")

                break  # valid input

            except ValueError as e:
                print(f"Invalid input: {e}")

        if self.graph is None:
            raise ValueError("Graph not created. Call create_graph() first.")
        
        # Create a copy of the graph to store influence zones
        influence_graph = self.graph.copy()
        
        # Initialize all nodes with no station assignment
        for node in influence_graph.nodes():
            influence_graph.nodes[node]['closest_station'] = None
            influence_graph.nodes[node][distance_col] = float('inf')
            influence_graph.nodes[node]['in_influence_zone'] = False
        
        # Initialize all edges with no station assignment
        for u, v in influence_graph.edges():
            influence_graph.edges[u, v]['closest_station'] = None
            influence_graph.edges[u, v]['cumulative_distance'] = 0.0
        
        # Get all station nodes
        station_nodes = [
            (node, data['station_id']) 
            for node, data in influence_graph.nodes(data=True) 
            if data.get('type') == 'station'
        ]
        
        print(f"Computing influence zones for {len(station_nodes)} stations...")
        print(f"Maximum upstream distance: {max_distance_km} km")
        
        max_distance_m = max_distance_km * 1000  # Convert to meters
        
        # Process each station
        for station_node, station_id in station_nodes:
            print(f"\nProcessing station {station_id}...")
            
            # Mark the station node itself
            influence_graph.nodes[station_node]['closest_station'] = station_id
            influence_graph.nodes[station_node][distance_col] = 0.0
            influence_graph.nodes[station_node]['in_influence_zone'] = True
            
            # Priority queue: (cumulative_distance, current_node, path_strahler_orders)
            # Using a list and sorting to prioritize lower Strahler orders
            queue = [(0.0, station_node, [])]
            visited = {station_node}
            
            while queue:
                # Sort by distance, then by max Strahler order (lower is better)
                queue.sort(key=lambda x: (x[0], max(x[2]) if x[2] else 0))
                current_dist, current_node, path_strahlers = queue.pop(0)
                
                # Get all upstream edges (predecessors)
                upstream_edges = [
                    (pred, current_node, influence_graph.edges[pred, current_node])
                    for pred in influence_graph.predecessors(current_node)
                ]
                
                # Sort by Strahler order (lower first) for prioritization
                upstream_edges.sort(key=lambda x: x[2]['strahler'])
                
                for pred_node, curr_node, edge_data in upstream_edges:
                    if pred_node in visited:
                        continue
                    
                    # Calculate cumulative distance
                    edge_length = edge_data['length']
                    new_dist = current_dist + edge_length
                    
                    # Check if within distance limit
                    if new_dist > max_distance_m:
                        continue
                    
                    # Check if this station is closer than any previously assigned
                    current_closest_dist = influence_graph.nodes[pred_node][distance_col]
                    
                    if new_dist < current_closest_dist:
                        # Update node attributes
                        influence_graph.nodes[pred_node]['closest_station'] = station_id
                        influence_graph.nodes[pred_node][distance_col] = new_dist
                        influence_graph.nodes[pred_node]['in_influence_zone'] = True
                        
                        # Update edge attributes
                        influence_graph.edges[pred_node, curr_node]['closest_station'] = station_id
                        influence_graph.edges[pred_node, curr_node]['cumulative_distance'] = new_dist
                        
                        # Mark as visited and add to queue
                        visited.add(pred_node)
                        
                        # Track Strahler orders in path
                        new_path_strahlers = path_strahlers + [edge_data['strahler']]
                        queue.append((new_dist, pred_node, new_path_strahlers))
            
            nodes_in_zone = sum(1 for n in visited)
            print(f"  Station {station_id}: {nodes_in_zone} nodes within {max_distance_km} km")
        
        # Store the influence graph as class attribute
        self.influence_graph = influence_graph
        
        # Print summary statistics
        print("\n=== Influence Zone Summary ===")
        assigned_nodes = sum(1 for n, d in influence_graph.nodes(data=True) 
                            if d['closest_station'] is not None)
        total_nodes = influence_graph.number_of_nodes()
        print(f"Nodes assigned to stations: {assigned_nodes}/{total_nodes} ({assigned_nodes/total_nodes*100:.1f}%)")
        
        return influence_graph


    def get_station_influence_geodataframe(self, station_id=None):
        """
        Get the influence zone for a station as a GeoDataFrame.
        
        Parameters:
        -----------
        station_id : str, optional
            Specific station ID. If None, returns all influence zones.
            
        Returns:
        --------
        gpd.GeoDataFrame
            Edges within the station's influence zone
        """
        if not hasattr(self, 'influence_graph') or self.influence_graph is None:
            raise ValueError("Influence zones not computed. Call compute_station_influence_zones() first.")
        
        edges_data = []
        
        for u, v, data in self.influence_graph.edges(data=True):
            closest_station = data.get('closest_station')
            
            # Filter by station_id if provided
            if station_id is not None and closest_station != station_id:
                continue
            
            if closest_station is not None:
                edge_dict = {
                    'geometry': data['geometry'],
                    'from_node': str(u),
                    'to_node': str(v),
                    'segment_id': data['segment_id'],
                    'strahler': data['strahler'],
                    'length': data['length'],
                    'closest_station': closest_station,
                    'cumulative_distance': data.get('cumulative_distance', 0.0),
                    'distance_km': data.get('cumulative_distance', 0.0) / 1000
                }
                edges_data.append(edge_dict)
        
        edges_gdf = gpd.GeoDataFrame(edges_data, crs=self.streams_gdf.crs)
        
        if station_id:
            print(f"Influence zone for station {station_id}: {len(edges_gdf)} edges")
        else:
            print(f"All influence zones: {len(edges_gdf)} edges")
        
        return edges_gdf