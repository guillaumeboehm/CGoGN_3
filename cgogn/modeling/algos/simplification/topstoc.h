/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: http://cgogn.unistra.fr/                                           *
* Contact information: cgogn@unistra.fr                                        *
*                                                                              *
*******************************************************************************/

#ifndef CGOGN_MODELING_ALGOS_TOPSTOC_H_
#define CGOGN_MODELING_ALGOS_TOPSTOC_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/io/surface/surface_import.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;

template <typename MESH>
void topstoc_vertex_selection(MESH &m, CellCache<MESH> &rg_chache, uint32 n)
{
	uint32 s = nb_cells<typename cgogn::mesh_traits<MESH>::Vertex>(m);
	uint32 i = 0;

	foreach_cell(m, [&](typename cgogn::mesh_traits<MESH>::Vertex v) -> bool { 
		if(i++ < s*0.95){
			rg_chache.template add<typename cgogn::mesh_traits<MESH>::Vertex>(v);
		}
		return true;
	});
}

template <typename MESH, typename Vertex, typename Face>
void compute_surface_data(const MESH &m, MESH &new_m,
						  const typename mesh_traits<MESH>::template Attribute<Vec3> *vertex_position,
						  const typename mesh_traits<MESH>::template Attribute<uint32> *vertex_anchor,
						  const CellMarkerStore<MESH, Vertex> &cm_selected,
						  const CellMarkerStore<MESH, Face> &cm_faces,
						  cgogn::io::SurfaceImportData &sd)
{
	std::vector<uint32> vertices = cm_selected.marked_cells();
	const uint32 nb_vertices = vertices.size();
	std::vector<uint32> faces = cm_faces.marked_cells();
	const uint32 nb_faces = faces.size();

	sd.reserve(nb_vertices, nb_faces);
	auto position = add_attribute<geometry::Vec3, CMap2::Vertex>(new_m, "position");

	// for storing the link between the prev vertices and the new ones
	std::unordered_map<uint32, uint32> v_map;

	//foreach cell of the prev mesh
	foreach_cell(m, [&](Vertex v) -> bool {
		if(cm_selected.is_marked(v)){ // if the vertex is to keep
			uint32 id = new_index<typename cgogn::mesh_traits<MESH>::Vertex>(new_m);
			(*position)[id] = value<Vec3>(m, vertex_position, v);
			std::cerr << "from:" << m.index_of(v) << " to:" << id;
			// if (value<uint32>(m, vertex_anchor, v) != m.index_of(v))
				std::cerr << " anchor: " << value<uint32>(m, vertex_anchor, v);
			std::cerr << std::endl;
			v_map.emplace(m.index_of(v), id);
			sd.vertices_id_.push_back(id);
		}
		return true;
	});

	//foreach face of the prev mesh
	foreach_cell(m, [&](Face f) -> bool {
		if(cm_faces.is_marked(f)){// if the face is to keep
			sd.faces_nb_vertices_.push_back(3);
			std::vector<typename cgogn::mesh_traits<MESH>::Vertex> iv = incident_vertices(m, f);
			sd.faces_vertex_indices_.push_back(v_map[value<uint32>(m, vertex_anchor, iv[0])]);
			sd.faces_vertex_indices_.push_back(v_map[value<uint32>(m, vertex_anchor, iv[1])]);
			sd.faces_vertex_indices_.push_back(v_map[value<uint32>(m, vertex_anchor, iv[2])]);

			// if (value<uint32>(m, vertex_anchor, iv[0]) != m.index_of(iv[0]) ||
			// 	value<uint32>(m, vertex_anchor, iv[1]) != m.index_of(iv[1]) ||
			// 	value<uint32>(m, vertex_anchor, iv[2]) != m.index_of(iv[2]))
			// 	std::cerr << "origin: " << m.index_of(iv[0]) << " "
			// 			  << m.index_of(iv[1]) << " "
			// 			  << m.index_of(iv[2])
			// 			  << "   new: " << v_map[value<uint32>(m, vertex_anchor, iv[0])] << " "
			// 			  << v_map[value<uint32>(m, vertex_anchor, iv[1])] << " "
			// 			  << v_map[value<uint32>(m, vertex_anchor, iv[2])] << std::endl;
		}
		return true;
	});
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
void topstoc(ui::MeshProvider<MESH> *mp, MESH &m, typename mesh_traits<MESH>::template Attribute<Vec3> *vertex_position, uint32 nb_vertices_to_keep)
{
	using Vertex = typename cgogn::mesh_traits<MESH>::Vertex;
	using Edge = typename cgogn::mesh_traits<MESH>::Edge;
	using Face = typename cgogn::mesh_traits<MESH>::Face;

	CellCache<MESH> rg_cache(m); //region growth cache


	//add vertex anchor attribute
	auto vertex_anchor = add_attribute<uint32, Vertex>(m, "anchor");

	//************ selection of vertices to keep
	topstoc_vertex_selection(m, rg_cache, nb_vertices_to_keep);

	//mark vertices with anchor selected vertex
	CellMarker<MESH, Vertex> cm_done(m);
	CellMarkerStore<MESH, Vertex> cm_selected(m);
	CellMarkerStore<MESH, Vertex> cm_rg_todo(m);
	CellMarkerStore<MESH, Face> cm_faces_to_keep(m);

	//mark the first selected vertices in the rg_cache 
	foreach_cell(rg_cache, [&](Vertex v) -> bool {
		cm_selected.mark(v);
		cm_rg_todo.mark(v);
		//mark its own vertex anchor if the vertex is selected
		value<uint32>(m, vertex_anchor, v) = m.index_of(v);
		return true;
	});


	//************* region-growth
	bool done = false;
	while(!done){
		done = true;
		foreach_cell(m, [&](Vertex v) -> bool {
			//if cell is not marked
			if(cm_rg_todo.is_marked(v) && !cm_done.is_marked(v))
			{
				done = false;
				std::cout << "ON VERTEX : " << m.index_of(v) << std::endl;
				//for each vertex in the one-ring
				foreach_adjacent_vertex_through_edge(m, v, [&](Vertex v2) -> bool {
					if(!cm_done.is_marked(v2) && !cm_selected.is_marked(v2))
					{
						//add to the rg_cache if the vertex is not marked
						cm_rg_todo.mark(v2);
						value<uint32>(m, vertex_anchor, v2) = value<uint32>(m, vertex_anchor, v);
						std::cout << "anchor : Vertex:" << m.index_of(v2) << " Parent:" << m.index_of(v) << " anchor:" << value<uint32>(m, vertex_anchor, v2) << std::endl;
					}
					else{
						std::cout << "pass : Vertex:" << m.index_of(v2) << " anchor:" << value<uint32>(m, vertex_anchor, v2) << std::endl;
					}
					return true;
				});
				//mark cell and anchor it
				cm_done.mark(v);
			}
			return true;
		});
	}




	////DEBUG
	foreach_cell(m, [&](Vertex v) -> bool {
		std::cout << "Vertex : " << m.index_of(v);
		if(!cm_selected.is_marked(v))
			std::cout << " not to keep";

		if (m.index_of(v) != value<uint32>(m, vertex_anchor, v))
			std::cout << " v:" << m.index_of(v) << " anchor:" << value<uint32>(m, vertex_anchor, v);
		std::cout << std::endl;
	});
	////DEBUG





	//store the faces to keep
	foreach_cell(m, [&](Face f) -> bool {
		std::vector<Vertex> iv = incident_vertices(m, f);

		if(value<uint32>(m, vertex_anchor, iv[0]) != value<uint32>(m, vertex_anchor, iv[1])
		&& value<uint32>(m, vertex_anchor, iv[0]) != value<uint32>(m, vertex_anchor, iv[2])
		&& value<uint32>(m, vertex_anchor, iv[1]) != value<uint32>(m, vertex_anchor, iv[2]))
			cm_faces_to_keep.mark(f);
		return true;
	});

	//************ simplification
	//create the new mesh that will be added to the mesh provider
	std::string name = "simplified_";
	mp->foreach_mesh([&](MESH *_m, const std::string& _name) -> bool {
		if (_m == &m)
			name += _name;
		return true;
	});

	MESH* new_m = mp->add_mesh(name);

	//setup surface_data
	cgogn::io::SurfaceImportData surface_data;
	compute_surface_data(m, *new_m, vertex_position, vertex_anchor.get(), cm_selected, cm_faces_to_keep, surface_data);

	//import surface data
	io::import_surface_data(*new_m, surface_data);

	std::shared_ptr<typename cgogn::mesh_traits<MESH>::template Attribute<Vec3>> new_vertex_position = cgogn::get_attribute<Vec3, Vertex>(*new_m, "position");

	if (new_vertex_position)
		mp->set_mesh_bb_vertex_position(new_m, new_vertex_position);

	// mp->emit_connectivity_changed(new_m);
	// mp->emit_attribute_changed(new_m, new_vertex_position.get());

	remove_attribute<Edge>(m, vertex_anchor);
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_TOPSTOC_H_
