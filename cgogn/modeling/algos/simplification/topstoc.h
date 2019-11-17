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
#include <cgogn/modeling/algos/subdivision.h>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;

template <typename MESH>
void topstoc_vertex_selection(MESH &m, CellCache<MESH> &selected, CellCache<MESH> &rg_chache, uint32 n)
{
	uint32 s = nb_cells(m);
	uint32 i = 0;
	uint32 step = s/n;
	foreach_cell(m, [&] (Vertex v) -> bool { 
		if(i%step == 0){
			selected.template add<Vertex>(v);
			rg_chache.template add<Vertex>(v);
		}
	})
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
void topstoc(MESH &m, typename mesh_traits<MESH>::template Attribute<Vec3> *vertex_position, uint32 nb_vertices_to_keep)
{
	using Vertex = typename cgogn::mesh_traits<MESH>::Vertex;
	using Edge = typename cgogn::mesh_traits<MESH>::Edge;
	using Face = typename cgogn::mesh_traits<MESH>::Face;

	CellCache<MESH> rg_cache(m); //region growth cache

	//add vertex anchor attribute
	auto vertex_anchor = add_attribute<uint32, Vertex>(m, "anchor");

	//selection of vertices to keep
	topstoc_vertex_selection(m, selected, rg_cache, nb_vertices_to_keep);

	//mark vertices with anchor selected vertex
	CellMarker<MESH, Vertex> cm_done(m);
	CellMarker<MESH, Vertex> cm_selected(m);

	//mark the first selected vertices in the rg_cache 
	foreach_cell(rg_cache, [&](Vertex v) -> bool {
		cm_selected.mark(v);
		return true;
	});

	
	//region-growth
	foreach_cell(rg_cache, [&](Vertex v) -> bool {
		//if cell is not marked
		if(!cm_done.is_marked(v))
		{
			//for each vertex in the one-ring
			foreach_incident_vertex(m, v, [&] (Vertex v2) -> bool {
				if(!cm.is_marked(v2) && !cm_selected.is_marked(v2)){
					//add to the rg_cache if the vertex is not marked
					rg_chache.template add<Vertex>(v2);
				}
				return true;
			});
			//mark cell and anchor it
			cm.mark(v);
			value<uint32>(m, vertex_anchor, v) = index_of(m, v);
		}
		return true;
	});

	//simplification

	// CellMarker<MESH, Face> cm_faces_to_keep(m);

	// foreach_cell(m, [&](Face f) -> bool {
	// 	vertices = incident_vertices(m, f);

	// 	// if the three vertices of the triangle have a different anchor
	// 	if(value<uint32>(m, vertex_anchor, vertices[0]) == value<uint32>(m, vertex_anchor, vertices[1]))
	// 		if(value<uint32>(m, vertex_anchor, vertices[1]) == value<uint32>(m, vertex_anchor, vertices[2]))
	// 			if(value<uint32>(m, vertex_anchor, vertices[0]) == value<uint32>(m, vertex_anchor, vertices[2])){
					
	// 				value<Vec3>(m, vertex_position, vertices[0]) = value<Vec3>(m, vertex_position, value<uint32>(m, vertex_anchor, vertices[0]));
	// 				value<Vec3>(m, vertex_position, vertices[1]) = value<Vec3>(m, vertex_position, value<uint32>(m, vertex_anchor, vertices[1]));
	// 				value<Vec3>(m, vertex_position, vertices[2]) = value<Vec3>(m, vertex_position, value<uint32>(m, vertex_anchor, vertices[2]));
	// 			}
		

	// 	return true;
	// });

}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_TOPSTOC_H_
