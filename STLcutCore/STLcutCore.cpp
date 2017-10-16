// STLcutCore.cpp : 定義 DLL 應用程式的匯出函式。
//

#include "stdafx.h"

#include <regex>
#include <xstring>
#include <windows.h>

void LoadStlBinary(const string& fileName, size_t& trianglecount, float*& dest)
{
	using namespace std;
	
	try
	{
		ifstream fin(fileName, ios::in | ios::binary);

		if (!fin.good())
		{
			throw new exception("檔名或路徑錯誤! 無法開啟檔案!");
		}

		fin.seekg(0, ios::end);
		string contents((unsigned int)fin.tellg(), '0');
		fin.seekg(0, ios::beg);
		fin.read(&contents[0], contents.size());
		fin.close();

		memcpy(&trianglecount, &contents[80], 4);

		char solid[6] = { '\0' };
		memcpy(solid, &contents[0], 5);
		if (strcmp(solid, "solid") == 0 && contents.size() != trianglecount * 50 + 84)
		{
			throw new domain_error("asci file.");
		}
		if (trianglecount < 0 || trianglecount > 1e8)
		{
			throw new out_of_range("triangle size incorrect.");
		}

		if (dest)
		{
			delete[] dest;
		}
		dest = new float[trianglecount * 3 * 3];

		float* tempp = dest;
		char* tempb = &contents[84] + 12;


		auto check_hash = [](float* f)
		{
			auto hashfunc = [](float* f)
			{
				return _Hash_seq((unsigned char*)f, 12);
			};
			auto hh1 = hashfunc(f);
			auto hh2 = hashfunc(f + 3);
			auto hh3 = hashfunc(f + 6);
			return hh1 == hh2 || hh2 == hh3 || hh3 == hh1;
		};
		auto check_inf_nan = [](float* f)
		{
			for (unsigned int j = 0; j < 9; j++)
			{
				if (isinf(*(f + j)) || isnan(*(f + j)))
				{
					return true;
				}
			}
			return false;
		};

		for (unsigned int i = 0; i < trianglecount; i++, tempb += 50)
		{
			memcpy(tempp, tempb, 36);

			if (check_inf_nan(tempp) || check_hash(tempp))
			{
				continue;
			}

			tempp += 9;
		}
		trianglecount = (tempp - dest) / 9;
	}
	catch (exception*)
	{
		trianglecount = 0;
		if (dest)
		{
			delete[] dest;
			dest = nullptr;
		}
	}
	return;
}

