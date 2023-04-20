#pragma once
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include "Particle.hpp"
#include "raylib-cpp.hpp"

#include <utility>

struct pair_hash
{
	template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2>& pair) const {
		return (27985032798461 * std::hash<T1>()(pair.first)) ^ std::hash<T2>()(pair.second);
	}
};

class Cell {
public:
	std::pair<int, int> key;
	raylib::Vector2 TR;
	raylib::Vector2 size;
	std::unordered_set<Particle*> particles;
	std::vector<Cell*> neighbours;

	Cell(float x, float y, float cell_size) {
		TR = raylib::Vector2(x, y);
		size = raylib::Vector2(cell_size, cell_size);
	}

	void RemoveParticle(Particle* p) {
		particles.erase(p);
	}

	void Draw() {
		DrawRectangleV(TR, size, GREEN);
	}

	void DrawNeighbours() {
		for (Cell* cell : neighbours) {
			cell->Draw();
		}
	}
};

class SpatialHash {
private:
	float cell_size;
	float width;
	float height;
	std::unordered_map<std::pair<int, int>, Cell*, pair_hash> table;

public:
	int numX, numY = 0;

	SpatialHash(float _cell_size, float _width, float _height): cell_size(_cell_size), width(_width), height(_height) {
		std::vector<Cell*> cells;
		numX = 0;
		int i, j;
		i = -1;
		for (float x = -cell_size; x < width+ _cell_size; x += _cell_size) {
			numY = 0;
			j = -1;
			for (float y = -cell_size; y < height+ _cell_size; y += _cell_size) {
				Cell* cell = new Cell(x, y, cell_size);
				cell->key = std::pair<int, int>(i,j);
				table.insert_or_assign(cell->key, cell);
				cells.push_back(cell);
				numY++;
				j++;
			}
			numX++;
			i++;
		}

		for (Cell* cell : cells) {
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					if (j != 0 || i != 0) {
						if (cell->key.first + i >= -1 && cell->key.second + j >= -1 && cell->key.first + i < floor(width / cell_size)+1 && cell->key.second + j < floor(height / cell_size)+1) {
							std::pair<int, int> key(cell->key.first + i, cell->key.second + j);
							cell->neighbours.push_back(table.at(key));
						}
					}
				}
			}
		}
	}

	Cell* GetCell(std::pair<int, int> k) {
		return table.at(k);
	}

	Cell* GetCell(Particle* p) {
		return table.at(p->cellKey);
	}

	void AddParticleToCell(Particle* p) {
		if (p->p.x < 0 || p->p.y < 0|| p->p.x > width || p->p.y > height || isnan(p->p.x-p->p.y)) return;
		if (p->cellKey != std::pair<int, int>(-1, -1)) table.at(p->cellKey)->RemoveParticle(p);
		std::pair<int, int> key = std::pair<int, int>((int)std::floor(p->p.x / cell_size), (int)std::floor(p->p.y / cell_size));
		Cell* cell = table.at(key);
		cell->particles.insert(p);
		p->cellKey = key;
	}
};