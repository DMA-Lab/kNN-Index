
#ifndef LCTD_FLATTEN_HASH_MAP_H
#define LCTD_FLATTEN_HASH_MAP_H

#include <vector>

namespace koala {
	template<typename kType, typename vType>
	class _iterator {
	
	};
	
	template<typename V>
	class my_openadd_hashmap {
#define LoadfactorInverse 1.5
#define nullele (4294967295U)
	public:
		//构造函数：初始化自定义类
		my_openadd_hashmap() : hash_table(0, std::make_pair(nullele, V())), _size(0) {}
		
		_iterator<unsigned int, V> begin() {
			return _iterator<unsigned int, V>();
		}
		
		_iterator<unsigned int, V> end() {
			return _iterator<unsigned int, V>();
		}
		
		unsigned int iterator() {
			unsigned int p = 0;
			while (p < hash_table.size() && hash_table[p].first == nullele)++p;
			return p;
		}
		
		unsigned int &next(unsigned int &p) {
			++p;
			while (p < hash_table.size() && hash_table[p].first == nullele)++p;
			return p;
		}
		
		bool has_next(unsigned int p) {
			return p != hash_table.size();
		}
		
		inline V &operator[](unsigned int key) {
			int idx = find(key);
			if (idx == -1) {
				idx = insert(key, V());
			}
			return hash_table[static_cast<unsigned int>(idx)].second;
		}
		
		inline std::pair<unsigned int, V> &get_with_idx(unsigned int i) { return hash_table[i]; }
		
		
		inline int insert(unsigned int key, V v) {
			if ((_size + 1) * LoadfactorInverse > hash_table.size()) {
				expand();
			}
			int idx = put(key, v, hash_table);
			if (idx >= 0)_size++;
			return idx;
		}
		
		inline int insert(std::pair<unsigned int, V> kv) {
			return insert(kv.first, kv.second);
		}
		
		/*
		 * if found, return idx of inner data structure
		 * if not found, return -1
		 */
		inline int find(unsigned int key) {
			if (_size == 0)
				return -1;
			unsigned int idx = key % hash_table.size();
			unsigned int tmpi = (idx + hash_table.size() - 1) % hash_table.size();
			while (hash_table[idx].first != nullele && hash_table[idx].first != key && tmpi != idx) {
				++idx;
				idx = idx == hash_table.size() ? 0 : idx;
			}
			if (hash_table[idx].first == key) {
				return static_cast<int>(idx);
			}
			return -1;
		}	
		
		inline void clear() {
			hash_table.clear(), hash_table.shrink_to_fit();
			//shrink_to_fit(): deque、vector 或 string 退回不需要的内存空间
			_size = 0;
		}
		
		inline std::vector<unsigned int> &release_data() {
			return hash_table;
		}
		
		inline unsigned int size() { return _size; }
	
	private:
		static inline int put(unsigned int key, V &v, std::vector<std::pair<unsigned int, V>> &table) {
			unsigned int idx = key % table.size();
			//std::cout<<"ini idx: ("<<idx<<", "<<key<<", "<<v<<")"<<std::endl;
			while (table[idx].first != nullele && table[idx].first != key) {
				++idx;
				idx = idx == table.size() ? 0 : idx;
			}
			//key 已 valid 存在
			//std::cout<<"after while idx: "<<idx<<std::endl;
			if (table[idx].first == nullele) {
				table[idx] = std::make_pair(key, v);
				return static_cast<int>(idx);
				//static_cast<type>(x) 把 idx 的类型强制转化成type类型 
			}
			// 是时候进行插入了

			return -1;
		}
		
		static inline int put(std::pair<unsigned int, V> p, std::vector<std::pair<unsigned int, V>> &table) {
			return put(p.first, p.second, table);
		}
		
		inline void expand() {
			unsigned int _new_size = 0;
			if (hash_table.size() == 0)_new_size = 1;
			else _new_size = hash_table.size() * 2;
			std::vector<std::pair<unsigned int, V>> newtable(_new_size, std::make_pair(nullele, V()));
			for (auto &a:hash_table)if (a.first != nullele)put(a, newtable);
			swap(newtable, hash_table);
		}
		
		std::vector<std::pair<unsigned int, V>> hash_table;
		unsigned int _size;
	};
}
#endif //LCTD_FLATTEN_HASH_MAP_H
