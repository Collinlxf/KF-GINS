#!/usr/bin/env python3
# coding: utf-8
# Export lane mapping results to GeoJSON for QGIS

import numpy as np
import json
import glob
import os
from pathlib import Path

print("=" * 60)
print("MonoLaneMapping GeoJSON Exporter for QGIS")
print("=" * 60)

# Lane type definitions
LANE_TYPES = {
    1: 'White Dash',
    2: 'White Solid',
    3: 'Double White',
    4: 'Yellow Dash',
    5: 'Yellow Solid',
    6: 'Double Yellow',
    20: 'Road Curb',
    21: 'Road Boundary',
}

LANE_COLORS = {
    1: '#808080',  # Gray
    2: '#FFFFFF',  # White
    3: '#FFFFFF',  # White
    4: '#FFFF00',  # Yellow
    5: '#FFFF00',  # Yellow
    6: '#FFFF00',  # Yellow
    20: '#FF0000',  # Red
    21: '#FF0000',  # Red
}

def export_final_map_to_geojson(map_file, output_file):
    """Export the final accumulated map from map.npy to GeoJSON"""
    
    print(f"\n[INFO] Loading map data: {map_file}")
    map_data = np.load(map_file, allow_pickle=True).item()
    
    if 'lanes_in_map' not in map_data:
        print("[ERROR] Invalid map file format")
        return False
    
    lanes = map_data['lanes_in_map']
    print(f"[INFO] Found {len(lanes)} lanes in map")
    
    # Create GeoJSON structure
    geojson = {
        "type": "FeatureCollection",
        "crs": {
            "type": "name",
            "properties": {
                "name": "urn:ogc:def:crs:EPSG::32633"  # UTM Zone 33N, adjust if needed
            }
        },
        "features": []
    }
    
    # Add each lane as a LineString feature
    for lane_id, lane_info in lanes.items():
        xyz = lane_info['xyz']
        category = lane_info.get('category', 0)
        ctrl_pts = lane_info.get('ctrl_pts', None)
        
        # Create LineString coordinates (X, Y, Z)
        coordinates = [[float(x), float(y), float(z)] for x, y, z in xyz]
        
        feature = {
            "type": "Feature",
            "properties": {
                "lane_id": int(lane_id),
                "category": int(category),
                "lane_type": LANE_TYPES.get(category, 'Unknown'),
                "color": LANE_COLORS.get(category, '#0000FF'),
                "num_points": len(xyz),
                "num_ctrl_pts": len(ctrl_pts) if ctrl_pts is not None else 0,
                "length_m": calculate_length(xyz),
            },
            "geometry": {
                "type": "LineString",
                "coordinates": coordinates
            }
        }
        
        geojson["features"].append(feature)
    
    # Save to file
    with open(output_file, 'w') as f:
        json.dump(geojson, f, indent=2)
    
    print(f"[SUCCESS] GeoJSON exported: {output_file}")
    print(f"[INFO] Total features: {len(geojson['features'])}")
    
    return True

def export_frame_to_geojson(json_file, output_file):
    """Export a single frame from JSON to GeoJSON"""
    
    print(f"\n[INFO] Loading frame: {json_file}")
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    if 'lane_lines' not in data:
        print("[ERROR] Invalid JSON format")
        return False
    
    lanes = data['lane_lines']
    print(f"[INFO] Found {len(lanes)} lanes in frame")
    
    # Create GeoJSON structure
    geojson = {
        "type": "FeatureCollection",
        "crs": {
            "type": "name",
            "properties": {
                "name": "urn:ogc:def:crs:EPSG::32633"
            }
        },
        "features": []
    }
    
    # Add each lane
    for idx, lane in enumerate(lanes):
        xyz = np.array(lane['xyz'])
        category = lane.get('category', 0)
        
        coordinates = [[float(x), float(y), float(z)] for x, y, z in xyz]
        
        feature = {
            "type": "Feature",
            "properties": {
                "lane_id": idx,
                "category": int(category),
                "lane_type": LANE_TYPES.get(category, 'Unknown'),
                "color": LANE_COLORS.get(category, '#0000FF'),
                "num_points": len(xyz),
                "length_m": calculate_length(xyz),
            },
            "geometry": {
                "type": "LineString",
                "coordinates": coordinates
            }
        }
        
        geojson["features"].append(feature)
    
    # Save to file
    with open(output_file, 'w') as f:
        json.dump(geojson, f, indent=2)
    
    print(f"[SUCCESS] GeoJSON exported: {output_file}")
    print(f"[INFO] Total features: {len(geojson['features'])}")
    
    return True

def export_all_frames_to_geojson(results_dir, output_dir):
    """Export all frames to separate GeoJSON files"""
    
    json_files = sorted(glob.glob(os.path.join(results_dir, '*.json')))
    print(f"\n[INFO] Found {len(json_files)} frames to export")
    
    os.makedirs(output_dir, exist_ok=True)
    
    for idx, json_file in enumerate(json_files):
        timestamp = Path(json_file).stem
        output_file = os.path.join(output_dir, f'{timestamp}.geojson')
        
        if (idx + 1) % 20 == 0:
            print(f"[INFO] Processing {idx + 1}/{len(json_files)}...")
        
        export_frame_to_geojson(json_file, output_file)
    
    print(f"[SUCCESS] Exported all {len(json_files)} frames to {output_dir}")

def calculate_length(xyz):
    """Calculate the length of a lane in meters"""
    if len(xyz) < 2:
        return 0.0
    
    distances = np.sqrt(np.sum(np.diff(xyz, axis=0)**2, axis=1))
    return float(np.sum(distances))

def export_trajectory_to_geojson(results_dir, output_file):
    """Export vehicle trajectory to GeoJSON"""
    
    print(f"\n[INFO] Extracting trajectory from frames...")
    
    json_files = sorted(glob.glob(os.path.join(results_dir, '*.json')))
    
    trajectory = []
    for json_file in json_files:
        # Assume vehicle is at origin (0, 0, 0) in each frame
        # You can extract actual pose if available
        trajectory.append([0.0, 0.0, 0.0])
    
    geojson = {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "properties": {
                "name": "Vehicle Trajectory",
                "num_points": len(trajectory),
                "total_frames": len(json_files)
            },
            "geometry": {
                "type": "LineString",
                "coordinates": trajectory
            }
        }]
    }
    
    with open(output_file, 'w') as f:
        json.dump(geojson, f, indent=2)
    
    print(f"[SUCCESS] Trajectory exported: {output_file}")

def create_qgis_style_file(output_file):
    """Create a QGIS style file (.qml) for better visualization"""
    
    qml_content = '''<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.0.0">
  <renderer-v2 type="categorizedSymbol" attr="lane_type">
    <categories>
      <category render="true" symbol="0" value="White Solid" label="White Solid"/>
      <category render="true" symbol="1" value="White Dash" label="White Dash"/>
      <category render="true" symbol="2" value="Yellow Solid" label="Yellow Solid"/>
      <category render="true" symbol="3" value="Yellow Dash" label="Yellow Dash"/>
      <category render="true" symbol="4" value="Road Curb" label="Road Curb"/>
    </categories>
    <symbols>
      <symbol type="line" name="0">
        <layer class="SimpleLine">
          <prop k="line_color" v="255,255,255,255"/>
          <prop k="line_width" v="1.5"/>
        </layer>
      </symbol>
      <symbol type="line" name="1">
        <layer class="SimpleLine">
          <prop k="line_color" v="200,200,200,255"/>
          <prop k="line_width" v="1.0"/>
          <prop k="line_style" v="dash"/>
        </layer>
      </symbol>
      <symbol type="line" name="2">
        <layer class="SimpleLine">
          <prop k="line_color" v="255,255,0,255"/>
          <prop k="line_width" v="1.5"/>
        </layer>
      </symbol>
      <symbol type="line" name="3">
        <layer class="SimpleLine">
          <prop k="line_color" v="255,255,0,255"/>
          <prop k="line_width" v="1.0"/>
          <prop k="line_style" v="dash"/>
        </layer>
      </symbol>
      <symbol type="line" name="4">
        <layer class="SimpleLine">
          <prop k="line_color" v="255,0,0,255"/>
          <prop k="line_width" v="2.0"/>
        </layer>
      </symbol>
    </symbols>
  </renderer-v2>
</qgis>'''
    
    with open(output_file, 'w') as f:
        f.write(qml_content)
    
    print(f"[INFO] QGIS style file created: {output_file}")

def main():
    # Find output directories
    output_dirs = glob.glob('outputs/*/visualization/segment-*/')
    results_dirs = glob.glob('outputs/*/results/segment-*/')
    
    if len(output_dirs) == 0 and len(results_dirs) == 0:
        print("[ERROR] No results found, please run ./run_demo.sh first")
        exit(1)
    
    # Create export directory
    if len(output_dirs) > 0:
        base_dir = Path(output_dirs[0]).parent.parent
    else:
        base_dir = Path(results_dirs[0]).parent.parent
    
    export_dir = os.path.join(base_dir, 'geojson')
    os.makedirs(export_dir, exist_ok=True)
    
    print(f"\n[INFO] Export directory: {export_dir}")
    
    # Export options
    print("\n" + "=" * 60)
    print("Choose export option:")
    print("=" * 60)
    print("  1) Final accumulated map (Recommended)")
    print("  2) Single frame (specify frame number)")
    print("  3) All frames (may take time)")
    print("  4) All of the above")
    print("")
    
    choice = input("Enter choice [1-4]: ").strip()
    
    if choice == '1' or choice == '4':
        # Export final map
        if len(output_dirs) > 0:
            map_files = glob.glob(os.path.join(output_dirs[0], 'map.npy'))
            if map_files:
                output_file = os.path.join(export_dir, 'final_map.geojson')
                export_final_map_to_geojson(map_files[0], output_file)
                
                # Create style file
                style_file = os.path.join(export_dir, 'final_map.qml')
                create_qgis_style_file(style_file)
    
    if choice == '2':
        # Export single frame
        if len(results_dirs) > 0:
            json_files = sorted(glob.glob(os.path.join(results_dirs[0], '*.json')))
            print(f"\nTotal frames: {len(json_files)}")
            frame_num = input("Enter frame number (1-{}): ".format(len(json_files))).strip()
            
            try:
                frame_idx = int(frame_num) - 1
                if 0 <= frame_idx < len(json_files):
                    timestamp = Path(json_files[frame_idx]).stem
                    output_file = os.path.join(export_dir, f'frame_{frame_num}_{timestamp}.geojson')
                    export_frame_to_geojson(json_files[frame_idx], output_file)
                else:
                    print("[ERROR] Invalid frame number")
            except ValueError:
                print("[ERROR] Invalid input")
    
    if choice == '3' or choice == '4':
        # Export all frames
        if len(results_dirs) > 0:
            frames_dir = os.path.join(export_dir, 'frames')
            export_all_frames_to_geojson(results_dirs[0], frames_dir)
    
    # Print summary
    print("\n" + "=" * 60)
    print("Export Complete!")
    print("=" * 60)
    print(f"\nExported files location: {export_dir}")
    print("\nHow to open in QGIS:")
    print("  1. Open QGIS")
    print("  2. Layer → Add Layer → Add Vector Layer")
    print("  3. Browse to: {}".format(export_dir))
    print("  4. Select final_map.geojson")
    print("  5. (Optional) Load final_map.qml for styling")
    print("\nLayer Properties:")
    print("  - Coordinate System: Custom (local coordinates)")
    print("  - Geometry: LineString with Z values")
    print("  - Attributes: lane_type, color, length_m, etc.")
    print("\n" + "=" * 60)

if __name__ == '__main__':
    main()
