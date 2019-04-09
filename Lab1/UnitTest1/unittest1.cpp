#include "stdafx.h"
#include "CppUnitTest.h"

#include<cmath>

#include "../Lab1/Header.h"
#include "../Lab1/Source.cpp"
#include "../Lab1/Generator.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTests
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		
		TEST_METHOD(TestMethod1)
		{
			int n = 4;
			vector<double> a = { 0, 1, 1, 1 };
			vector<double> b = { -2, -2, -2, -2 };
			vector<double> c = { 1, 1, 1, 0 };
			vector<double> d = { 0.04, 0.04, 0.04, 0.04 };

			vector<double> sol = { -0.08, -0.12, -0.12, -0.08 };
			vector<double> u(n);

			thomas(a, b, c, u, d);

			for (int i = 0; i < n; i++) {
				Assert::AreEqual(u[i], sol[i]);
			}
		}

		TEST_METHOD(TestMethod2)
		{
			int n = 3;
			vector<double> a = { 0, 4, 3 };
			vector<double> b = { 9, -7, 8};
			vector<double> c = { 1, 2, 0 };
			vector<double> d = { 5, 6, 2 };

			vector<double> sol = { 0.6, -0.4, 0.4 };
			vector<double> u(n);

			thomas(a, b, c, u, d);

			for (int i = 0; i < n; i++) {
				Assert::IsTrue(std::abs(u[i] - sol[i]) <= 0.001);
			}
		}

		TEST_METHOD(TestMethod3)
		{
			int n = 4;
			vector<double> a = { 0, 3, 6, 9};
			vector<double> b = { 1, 4, 7, 10 };
			vector<double> c = { 2, 5, 8, 0 };
			vector<double> d = { 3, 12, 21, 19 };

			vector<double> sol = { 1, 1, 1, 1 };
			vector<double> u(n);

			thomas(a, b, c, u, d);

			for (int i = 0; i < n; i++) {
				Assert::IsTrue(std::abs(u[i] - sol[i]) <= 0.001);
			}
		}

		TEST_METHOD(Generated1)
		{
			int n = 4;
			vector<double> a(n);
			vector<double> b(n);
			vector<double> c(n);
			vector<double> d(n);
			vector<double> u(n);
			generate_thomas(n, a, b, c, u, d);
			thomas(a, b, c, u, d);

			for (int i = 0; i < n; i++) {
				Assert::IsTrue(std::abs(u[i] - 1) <= 0.001);
			}
		}

		TEST_METHOD(With10kElems)
		{
			int n = 10000;
			vector<double> a(n);
			vector<double> b(n);
			vector<double> c(n);
			vector<double> d(n);
			vector<double> u(n);
			generate_thomas(n, a, b, c, u, d);
			thomas(a, b, c, u, d);

			for (int i = 0; i < n; i++) {
				Assert::IsTrue(std::abs(u[i] - 1) <= 0.001);
			}
		}

		TEST_METHOD(With100kElems)
		{
			int n = 100000;
			vector<double> a(n);
			vector<double> b(n);
			vector<double> c(n);
			vector<double> d(n);
			vector<double> u(n);
			generate_thomas(n, a, b, c, u, d);
			thomas(a, b, c, u, d);

			for (int i = 0; i < n; i++) {
				Assert::IsTrue(std::abs(u[i] - 1) <= 0.001);
			}
		}

		TEST_METHOD(With1MilElems)
		{
			int n = 1000000;
			vector<double> a(n);
			vector<double> b(n);
			vector<double> c(n);
			vector<double> d(n);
			vector<double> u(n);
			generate_thomas(n, a, b, c, u, d);
			thomas(a, b, c, u, d);

			for (int i = 0; i < n; i++) {
				Assert::IsTrue(std::abs(u[i] - 1) <= 0.001);
			}
		}


		TEST_METHOD(CyclicRed)
		{

			int n = 7;
			vector<double> a(n);
			vector<double> b(n);
			vector<double> c(n);
			vector<double> d(n);
			vector<double> u(n);
			generate_thomas(n, a, b, c, u, d);
			cyclic_reduction_omp(a, b, c, u, d);

			for (int i = 0; i < n; i++) {
				Assert::IsTrue(std::abs(u[i] - 1) <= 0.001);
			}
		}
	};
}