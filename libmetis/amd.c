#include "metislib.h"

// helpers

void print_csr(idx_t nvtxs, idx_t* xadj, idx_t* adjncy)
{
	idx_t i, j;

	for (i = 0; i < nvtxs; i++)
	{
		printf("xadj[%d]: %d\n", i, xadj[i]);
	}
	for (i = 0; i < nvtxs; i++)
	{
		for (j = xadj[i]; j < xadj[i + 1]; j++)
		{
			printf("adjncy[%d]: %d\n", j, adjncy[j]);
		}
	}
}

void print_graph(my_graph_t* graph)
{
	for (idx_t i = 0; i < graph->n; i++)
	{
		printf("first[%d]: %d\n", i, graph->first[i]);
	}
	for (idx_t i = 0; i < graph->m; i++)
	{
		printf("nx[%d]: %d\n", i, graph->nx[i]);
	}
	for (idx_t i = 0; i < graph->m; i++)
	{
		printf("prev[%d]: %d\n", i, graph->prev[i]);
	}
	for (idx_t i = 0; i < graph->m; i++)
	{
		printf("adj[%d]: %d\n", i, graph->adj[i]);
	}
}

void print_node_adj(my_graph_t* graph, idx_t node)
{
	printf("node: %d\n", node);

	for (idx_t i = graph->first[node];i != -1;i = graph->nx[i])
	{
		printf("%d ", graph->adj[i]);
	}

	printf("\n");
}

void print_graph_adj(my_graph_t* graph, char* str)
{
	printf("\n %s \n", str);
	for (idx_t i = 0;i < graph->n;i++)
	{
		print_node_adj(graph, i);
	}
}

void print_mask(idx_t* mask, idx_t n, char* str)
{
	printf("\n %s \n", str);
	for (idx_t i = 0;i < n;i++)
	{
		printf("%d ", mask[i]);
	}
	printf("\n");
}

int is_valid(my_graph_t* graph)
{
	for (idx_t i = 0; i < graph->m; i++)
	{
		if (graph->nx[i] != -1 && graph->prev[graph->nx[i]] != i)
		{
			return 0;
		}
	}

	return 1;
}

// common functions

void add_edge(my_graph_t* graph, idx_t u, idx_t v)
{
	graph->adj[graph->m] = v;
	graph->nx[graph->m] = graph->first[u];
	graph->first[u] = graph->m;
	graph->inv_first[graph->m] = u;
	graph->prev[graph->m] = -1;
	if (graph->nx[graph->m] != -1)
	{
		graph->prev[graph->nx[graph->m]] = graph->m;
	}
	graph->counts[u] += 1;

	graph->m += 1;

	if (graph->m == graph->m_max)
	{
		printf("add_edge: memory not enough \n");
		exit(0);

		graph->m_max = 2 * graph->m;

		graph->nx = irealloc(graph->nx, graph->m_max, "add_edge: nx");
		graph->adj = irealloc(graph->adj, graph->m_max, "add_edge: adj");
		graph->prev = irealloc(graph->prev, graph->m_max, "add_edge: prev");
		graph->inv_first = irealloc(graph->inv_first, graph->m_max, "add_edge: inv_first");
	}
}

//удаление конкретного элемента в списке смежности
void remove_edge(my_graph_t* graph, idx_t u, idx_t v)
{
	idx_t index = graph->first[u];
	if (graph->adj[index] == v)
	{
		if (graph->nx[index] != -1)
		{
			graph->prev[graph->nx[index]] = graph->prev[graph->first[u]];
			graph->inv_first[graph->nx[index]] = u;
		}
		graph->first[u] = graph->nx[index];
	}
	else
	{
		idx_t prev_index = -1;
		while (graph->adj[index] != v)
		{
			prev_index = index;
			index = graph->nx[index];

			if (index == -1)
			{
				printf("remove_edge: adjacency node not found \n");
				exit(0);
			}
		}
		if (graph->nx[index] != -1)
		{
			graph->prev[graph->nx[index]] = graph->prev[graph->nx[prev_index]];
		}
		graph->nx[prev_index] = graph->nx[index];
	}

	idx_t last_index = graph->m - 1;
	if (index != last_index)
	{
		// заполняем освободившийся индекс adj
		graph->adj[index] = graph->adj[last_index];

		// заполняем освободившийся индекс nx
		if (graph->nx[last_index] != -1)
		{
			graph->prev[graph->nx[last_index]] = index;
		}
		graph->nx[index] = graph->nx[last_index];

		// меняем индекс который ссылается на перемещённый индекс
		if (graph->prev[last_index] != -1)
		{
			graph->prev[index] = graph->prev[graph->nx[graph->prev[last_index]]];
			graph->nx[graph->prev[last_index]] = index;
		}
		else
		{
			idx_t first_index = graph->inv_first[last_index];
			graph->prev[index] = graph->prev[graph->first[first_index]];
			graph->first[first_index] = index;
			graph->inv_first[index] = first_index;
		}
	}

	graph->counts[u] -= 1;
	graph->m -= 1;
}

// mmd algorithm

void mmd_init(idx_t* xadj, idx_t* adjncy, my_graph_t* graph, idx_t* degrees, idx_t* perm, idx_t* iperm)
{
	for (idx_t i = 0; i < graph->n; i++)
	{
		graph->first[i] = -1;
		graph->counts[i] = 0;

		perm[i] = -1;
		iperm[i] = -1;

		degrees[i] = xadj[i + 1] - xadj[i];

		for (idx_t j = xadj[i]; j < xadj[i + 1]; j++)
		{
			add_edge(graph, i, adjncy[j]);
			// не добавляем обратное ребро, посольку в metis уже симетричная матрица на входе
		}
	}

	//only for test
	/*add_edge(graph, 0, 3);add_edge(graph, 3, 0);
	add_edge(graph, 1, 2);add_edge(graph, 2, 1);
	add_edge(graph, 1, 4);add_edge(graph, 4, 1);
	add_edge(graph, 1, 5);add_edge(graph, 5, 1);
	add_edge(graph, 2, 3);add_edge(graph, 3, 2);
	add_edge(graph, 2, 4);add_edge(graph, 4, 2);
	add_edge(graph, 5, 0);add_edge(graph, 0, 5);
	degrees[0] = 2;
	degrees[1] = 3;
	degrees[2] = 3;
	degrees[3] = 2;
	degrees[4] = 2;
	degrees[5] = 2;*/
}

void eliminate_node(my_graph_t* graph, idx_t node, idx_t* degrees, idx_t* mask)
{
	for (idx_t i = graph->first[node];i != -1;i = graph->nx[i])
	{
		idx_t adjacency_node = graph->adj[i];

		// помечаем все смежные вершины вершине adjacency_node
		iset(graph->n, 0, mask);
		for (idx_t j = graph->first[adjacency_node];j != -1;j = graph->nx[j])
		{
			mask[graph->adj[j]] = 1;
		}

		// если не хватает рёбер для клики то добавляем их
		for (idx_t j = graph->first[node];j != -1;j = graph->nx[j])
		{
			if (mask[graph->adj[j]] != 1 && graph->adj[j] != adjacency_node)
			{
				add_edge(graph, adjacency_node, graph->adj[j]);
				add_edge(graph, graph->adj[j], adjacency_node);

				degrees[adjacency_node] += 1;
				degrees[graph->adj[j]] += 1;
			}
		}
	}

	// удаляем все рёбра инцидентные вершине node
	idx_t index = graph->first[node];
	while (index != -1)
	{
		idx_t adjacency_node = graph->adj[index];
		index = graph->nx[index];

		remove_edge(graph, node, adjacency_node);
		remove_edge(graph, adjacency_node, node);

		degrees[node] -= 1;
		degrees[adjacency_node] -= 1;
	}
}

void mmd(idx_t* xadj, idx_t* adjncy, idx_t* perm, idx_t* iperm, my_graph_t* graph, idx_t* degrees, idx_t* mask)
{
	//only for test
	//graph->n = 6;

	mmd_init(xadj, adjncy, graph, degrees, perm, iperm);

	idx_t index = 0;
	while (graph->n != index)
	{
		idx_t min_degree_node = -1;
		idx_t min_degree = -1;

		//todo use bucketsort, ikvsort
		for (idx_t i = 0;i < graph->n;i++)
		{
			// выбранные вершины пропускаем
			if (iperm[i] != -1)
			{
				continue;
			}

			if (min_degree < 0 || min_degree > degrees[i])
			{
				min_degree = degrees[i];
				min_degree_node = i;
			}
		}

		//print_graph_adj(graph, "graph:");

		eliminate_node(graph, min_degree_node, degrees, mask);

		perm[index] = min_degree_node;
		iperm[min_degree_node] = index;

		index++;
	}
}

// amd algorithm

void amd_init(idx_t* xadj, idx_t* adjncy, my_graph_t* graph, my_graph_t* graph_e, my_graph_t* graph_l, my_graph_t* graph_h, my_graph_t* graph_s, idx_t* degrees, idx_t* perm, idx_t* iperm, idx_t* mask_lp, idx_t* mask_ep)
{
	for (idx_t i = 0; i < graph->n; i++)
	{
		graph->first[i] = -1;
		graph_e->first[i] = -1;
		graph_l->first[i] = -1;
		graph_h->first[i] = -1;
		graph_s->first[i] = -1;

		graph->counts[i] = 0;
		graph_e->counts[i] = 0;
		graph_l->counts[i] = 0;
		graph_h->counts[i] = 0;
		graph_s->counts[i] = 0;

		perm[i] = -1;
		iperm[i] = -1;

		degrees[i] = xadj[i + 1] - xadj[i];

		mask_lp[i] = -1;
		mask_ep[i] = -1;

		add_edge(graph_s, i, i);

		for (idx_t j = xadj[i]; j < xadj[i + 1]; j++)
		{
			add_edge(graph, i, adjncy[j]);
			// не добавляем обратное ребро, посольку в metis уже симетричная матрица на входе
		}
	}

	//only for test
	/*add_edge(graph, 0, 3);add_edge(graph, 3, 0);
	add_edge(graph, 0, 5);add_edge(graph, 5, 0);
	add_edge(graph, 1, 4);add_edge(graph, 4, 1);
	add_edge(graph, 1, 5);add_edge(graph, 5, 1);
	add_edge(graph, 1, 8);add_edge(graph, 8, 1);
	add_edge(graph, 2, 4);add_edge(graph, 4, 2);
	add_edge(graph, 2, 5);add_edge(graph, 5, 2);
	add_edge(graph, 2, 6);add_edge(graph, 6, 2);
	add_edge(graph, 3, 6);add_edge(graph, 6, 3);
	add_edge(graph, 3, 7);add_edge(graph, 7, 3);
	add_edge(graph, 4, 6);add_edge(graph, 6, 4);
	add_edge(graph, 4, 8);add_edge(graph, 8, 4);
	add_edge(graph, 6, 7);add_edge(graph, 7, 6);
	add_edge(graph, 6, 8);add_edge(graph, 8, 6);
	add_edge(graph, 6, 9);add_edge(graph, 9, 6);
	add_edge(graph, 7, 8);add_edge(graph, 8, 7);
	add_edge(graph, 7, 9);add_edge(graph, 9, 7);
	add_edge(graph, 8, 9);add_edge(graph, 9, 8);*/

	/*degrees[0] = 2;
	degrees[1] = 3;
	degrees[2] = 3;
	degrees[3] = 3;
	degrees[4] = 4;
	degrees[5] = 3;
	degrees[6] = 6;
	degrees[7] = 4;
	degrees[8] = 5;
	degrees[9] = 3;*/

	/*degrees[0] = 0;
	degrees[1] = 1;
	degrees[2] = 2;
	degrees[3] = 3;
	degrees[4] = 4;
	degrees[5] = 5;
	degrees[6] = 6;
	degrees[7] = 7;
	degrees[8] = 8;
	degrees[9] = 9;*/
}

void fill_w(my_graph_t* graph_e, my_graph_t* graph_l, my_graph_t* graph_s, idx_t node, idx_t* w)
{
	for (idx_t i = graph_e->first[node];i != -1;i = graph_e->nx[i])
	{
		if (w[graph_e->adj[i]] < 0)
		{
			w[graph_e->adj[i]] = graph_l->counts[graph_e->adj[i]];
		}

		w[graph_e->adj[i]] -= graph_s->counts[node];
	}
}

idx_t is_indistinguishable(my_graph_t* graph, my_graph_t* graph_e, idx_t* mask, idx_t principal_node, idx_t none_principal_node)
{
	if (graph->counts[principal_node] != graph->counts[none_principal_node] || graph_e->counts[principal_node] != graph_e->counts[none_principal_node])
	{
		return 0;
	}

	/*print_node_adj(graph, principal_node);
	print_node_adj(graph, none_principal_node);*/

	for (idx_t i = graph->first[none_principal_node];i != -1;i = graph->nx[i])
	{
		if (mask[graph->adj[i]] != 1)
		{
			return 0;
		}
	}

	/*print_node_adj(graph_e, principal_node);
	print_node_adj(graph_e, none_principal_node);*/

	for (idx_t i = graph_e->first[none_principal_node];i != -1;i = graph_e->nx[i])
	{
		if (mask[graph_e->adj[i]] != 2)
		{
			return 0;
		}
	}

	return 1;
}

void amd_eliminate_node(my_graph_t* graph, my_graph_t* graph_e, my_graph_t* graph_l, my_graph_t* graph_h, my_graph_t* graph_s, idx_t node, idx_t* degrees, idx_t* mask_lp, idx_t* mask_ep, idx_t* mask, idx_t k, idx_t* w)
{
	idx_t lp_count = 0;
	iset(graph->n, -1, w);
	// находим Lp, Ep с помощью mask_lp, mask_ep
	for (idx_t i = graph->first[node];i != -1;i = graph->nx[i])
	{
		fill_w(graph_e, graph_l, graph_s, graph->adj[i], w);
		lp_count++;
		mask_lp[graph->adj[i]] = node;
	}
	for (idx_t i = graph_e->first[node];i != -1;i = graph_e->nx[i])
	{
		idx_t adjacency_element = graph_e->adj[i];
		mask_ep[adjacency_element] = node;

		for (idx_t j = graph_l->first[adjacency_element];j != -1;j = graph_l->nx[j])
		{
			// можем несколько раз помечать одну и ту же переменную
			if (mask_lp[graph_l->adj[j]] != node)
			{
				fill_w(graph_e, graph_l, graph_s, graph->adj[j], w);
				lp_count++;
			}
			mask_lp[graph_l->adj[j]] = node;
		}
	}
	mask_lp[node] = -1;

	/*print_mask(mask_lp, graph->n, "mask_lp:");
	print_mask(mask_ep, graph->n, "mask_ep:");*/

	for (idx_t i = 0; i < graph->n; i++)
	{
		if (mask_lp[i] == node)
		{
			// преобразовываем Ai
			idx_t index = graph->first[i];
			while (index != -1)
			{
				idx_t adjacency_node = graph->adj[index];
				index = graph->nx[index];

				if (mask_lp[adjacency_node] == node || adjacency_node == node)
				{
					remove_edge(graph, i, adjacency_node);
					remove_edge(graph, adjacency_node, i);
				}
			}

			idx_t l_sum = 0;
			// преобразовываем Ei, Li
			index = graph_e->first[i];
			while (index != -1)
			{
				idx_t adjacency_element = graph_e->adj[index];
				index = graph_e->nx[index];

				if (mask_ep[adjacency_element] == node)
				{
					remove_edge(graph_e, i, adjacency_element);
					remove_edge(graph_l, adjacency_element, i);
				}
				// учитываем только те элементы которые останутся после преобразования
				else
				{
					l_sum += w[adjacency_element] > 0 ? w[adjacency_element] : graph_l->counts[adjacency_element];
				}
			}
			add_edge(graph_e, i, node);
			add_edge(graph_l, node, i);

			// поскольку для суперпеременных мы уже удалили их списки смежности, то будет только одно лишнее вхождение i в lp_count
			// а также ни одного вхождения i в graph->counts[i]
			degrees[i] = min(min(graph->n - k, degrees[i] + lp_count - 1), graph->counts[i] + lp_count + l_sum - 1);

			add_edge(graph_h, (graph->counts[i] + graph_e->counts[i]) % (graph->n - 1) + 1, i);
		}
	}

	//print_graph_adj(graph_h, "graph_h:");

	// удаляем данные для найденных суперпеременных
	for (idx_t i = 0; i < graph_h->n; i++)
	{
		idx_t index_h = graph_h->first[i];
		if (index_h == -1)
		{
			continue;
		}

		idx_t principal_node = graph_h->adj[index_h];

		// заполняем маску смежных вершин в графе разбиений
		iset(graph->n, 0, mask);
		for (idx_t j = graph->first[principal_node];j != -1;j = graph->nx[j])
		{
			mask[graph->adj[j]] = 1;
		}
		for (idx_t j = graph_e->first[principal_node];j != -1;j = graph_e->nx[j])
		{
			if (mask[graph_e->adj[j]] == 1)
			{
				printf("amd_eliminate_node: graph intercepts graph_e \n");
				exit(0);
			}
			mask[graph_e->adj[j]] = 2;
		}

		index_h = graph_h->nx[index_h];
		while (index_h != -1)
		{
			idx_t none_principal_node = graph_h->adj[index_h];
			index_h = graph_h->nx[index_h];

			if (is_indistinguishable(graph, graph_e, mask, principal_node, none_principal_node))
			{
				degrees[principal_node] -= graph_s->counts[none_principal_node];

				// преобразовываем Aj, с j=none_principal_node
				idx_t index = graph->first[none_principal_node];
				while (index != -1)
				{
					idx_t adjacency_node = graph->adj[index];
					index = graph->nx[index];

					remove_edge(graph, none_principal_node, adjacency_node);
					remove_edge(graph, adjacency_node, none_principal_node);
				}

				// преобразовываем Ej, с j=none_principal_node
				index = graph_e->first[none_principal_node];
				while (index != -1)
				{
					idx_t adjacency_element = graph_e->adj[index];
					index = graph_e->nx[index];

					remove_edge(graph_e, none_principal_node, adjacency_element);
					remove_edge(graph_l, adjacency_element, none_principal_node);
				}

				// преобразовываем списки суперпеременных
				index = graph_s->first[none_principal_node];
				while (index != -1)
				{
					idx_t adjacency_node = graph_s->adj[index];
					index = graph_s->nx[index];

					remove_edge(graph_s, none_principal_node, adjacency_node);
					add_edge(graph_s, principal_node, adjacency_node);
				}
			}

			remove_edge(graph_h, i, none_principal_node);
		}
		remove_edge(graph_h, i, principal_node);
	}

	// очищаем Ep, Ap к этому моменту должно быть уже очищено
	idx_t index = graph_e->first[node];
	while (index != -1)
	{
		idx_t adjacency_element = graph_e->adj[index];
		index = graph_e->nx[index];

		remove_edge(graph_e, node, adjacency_element);
		remove_edge(graph_l, adjacency_element, node);
	}
}

void amd(idx_t* xadj, idx_t* adjncy, idx_t* perm, idx_t* iperm, my_graph_t* graph, my_graph_t* graph_e, my_graph_t* graph_l, my_graph_t* graph_h, my_graph_t* graph_s, idx_t* degrees, idx_t* mask_lp, idx_t* mask_ep, idx_t* mask, idx_t* w)
{
	//only for test
	/*graph->n = 10;
	graph_e->n = 10;
	graph_l->n = 10;
	graph_h->n = 10;
	graph_s->n = 10;*/

	amd_init(xadj, adjncy, graph, graph_e, graph_l, graph_h, graph_s, degrees, perm, iperm, mask_lp, mask_ep);

	idx_t index = 0;
	while (graph->n != index)
	{
		idx_t min_degree_node = -1;
		idx_t min_degree = -1;

		//todo use bucketsort, ikvsort
		for (idx_t i = 0;i < graph->n;i++)
		{
			// выбранные вершины или вершины в суперпеременных пропускаем
			if (iperm[i] != -1 || graph_s->counts[i] == 0)
			{
				continue;
			}

			if (min_degree < 0 || min_degree > degrees[i])
			{
				min_degree = degrees[i];
				min_degree_node = i;
			}
		}

		amd_eliminate_node(graph, graph_e, graph_l, graph_h, graph_s, min_degree_node, degrees, mask_lp, mask_ep, mask, index, w);

		/*print_graph_adj(graph, "graph:");
		print_graph_adj(graph_e, "graph_e:");
		print_graph_adj(graph_l, "graph_l:");*/

		//print_graph_adj(graph_s, "graph_s:");

		for (idx_t i = graph_s->first[min_degree_node];i != -1;i = graph_s->nx[i])
		{
			perm[index] = graph_s->adj[i];
			iperm[graph_s->adj[i]] = index;
			index++;
		}
	}
}