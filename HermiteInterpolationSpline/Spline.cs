using System;
using System.Diagnostics;

namespace HermiteInterpolationSpline
{
   internal class Spline
   {
      private int SourcePointNumber;
      private double[]? _sourceX;
      private double[]? _sourceF;
      private double[]? _sourceDF;
      private int[]? _subareaSplits;

      private List<double> _meshX;
      private List<double> _meshSpline;

      public Spline() 
      {
         _meshX = new List<double>();
         _meshSpline = new List<double>();
      }

      public void InputSourcePoints(string pointsPath, string parametersPath)
      {
         try
         {
            using (StreamReader sr = new StreamReader(pointsPath))
            {
               _sourceX = sr.ReadLine().Split().Select(double.Parse).ToArray();
               _sourceF = sr.ReadLine().Split().Select(double.Parse).ToArray();
               SourcePointNumber = _sourceF.Length;
               _sourceDF = new double[SourcePointNumber];
            }

            if (_sourceX == null || _sourceX.Length == 0)
            {
               throw new ArgumentException("Ошибка при чтении координатых точек.");
            }

            if (_sourceX.Length < 3)
            {
               throw new ArgumentException("Недостаточно исходных данных. Для построения сплайна требуется 3 точки. ");
            }

            if (_sourceF == null || _sourceX.Length != _sourceF.Length) 
            { 
               throw new ArgumentException("Ошибка при чтении значений функции в точках."); 
            }

            SourcePointNumber = _sourceX.Length;

            using (StreamReader sr = new StreamReader(parametersPath))
            {
               _subareaSplits = sr.ReadLine().Split().Select(int.Parse).ToArray();
            }

            if (_subareaSplits == null || _subareaSplits.Length < _sourceX.Length - 1)
            {
               throw new ArgumentException("Ошибка при чтении количества разбиений подобластей сплайна.");
            }

            if (_subareaSplits.Length < 1)
            {
               throw new ArgumentException("Количество разбиений области должно быть больше или равным единице.");
            }
         }
         catch(Exception ex) 
         {
            Console.WriteLine(ex.Message);
            Process.GetCurrentProcess().Kill();
         }

      }

      public void Compute()
      {
         SetMesh();
         Derivation();

         List<int> indexes = new() { 0 };
         for (int i = 0; i < SourcePointNumber - 1; i++)
            indexes.Add(indexes[i] + _subareaSplits[i]);

         double Xi = 0;
         double hi = 0;
         _meshSpline.Add(_sourceF[0]);
         for (int i = 0; i < SourcePointNumber - 1; i++)
         {
            hi = _sourceX[i + 1] - _sourceX[i];
            for (int j = indexes[i]; j < indexes[i + 1]; j++)
            {
               Xi = ( _meshX[j + 1] - _meshX[ indexes[i] ] ) / hi;
               _meshSpline.Add
                  (
                  _sourceF[i] * HermiteBasis.Psi0(Xi) 
                  + _sourceF[i + 1] * HermiteBasis.Psi2(Xi)
                  + _sourceDF[i] * HermiteBasis.Psi1(Xi, hi) 
                  + _sourceDF[i + 1] * HermiteBasis.Psi3(Xi, hi)
                  );
            }
         }
      }

      public void Output(string outputPath)
      {
         using (StreamWriter sw = new StreamWriter(outputPath))
         {
            for (int i = 0; i < _meshX.Count; i++)
               sw.Write("{0:e15} ", _meshX[i]);

            sw.WriteLine();

            for (int i = 0; i < _meshSpline.Count; i++)
               sw.Write("{0:e15} ", _meshSpline[i]);
         }
      }

      // Создание координатной сетки.
      private void SetMesh()
      {
         double step;

         for (int i = 0; i < SourcePointNumber - 1; i++)
         {
            step = (_sourceX[i + 1] - _sourceX[i]) / _subareaSplits[i];
            for (int j = 0; j < _subareaSplits[i]; j++)
               _meshX.Add(_sourceX[i] + j * step);
         }
         _meshX.Add(_sourceX[SourcePointNumber - 1]);
      }

      // Аппроксимация производных исходной функции по трём точкам
      // с использованием квадратичного полинома Лагранжа.
      private void Derivation()
      {
         int n = _sourceF.Length;
         double hPrev;
         double hNext;

         // Левая производная.
         hPrev = _sourceX[1] - _sourceX[0];
         hNext = _sourceX[2] - _sourceX[1];
         _sourceDF[0] 
            = _sourceF[0] * -(2 * hPrev + hNext) / hPrev / (hPrev + hNext)
            + _sourceF[1] * (hPrev + hNext) / hPrev / hNext 
            + _sourceF[2] * -hPrev / hNext / (hPrev + hNext);

         // Центральная производная.
         hPrev = _sourceX[1] - _sourceX[0];
         for (int i = 1; i < n - 1; i++)
         {
            hNext = _sourceX[i + 1] - _sourceX[i];
            _sourceDF[i] 
               = _sourceF[i - 1] * Dl(hPrev, hNext)
               + _sourceF[i] * Dc(hPrev, hNext)
               + _sourceF[i + 1] * Dr(hPrev, hNext);

            hPrev = hNext;
         }

         // Правая производная.
         hPrev = _sourceX[n - 2] - _sourceX[n - 3];
         hNext = _sourceX[n - 1] - _sourceX[n - 2];
         _sourceDF[n - 1] =
            _sourceF[n - 3] * hPrev / hPrev / (hPrev + hNext)
            + _sourceF[n - 2] * -(hPrev + hNext) / hPrev / hNext
            + _sourceF[n - 1] * (2 * hNext + hPrev) / hNext / (hPrev + hNext);
      }

      private static double Dl(double hPrev, double hNext)
         => -hNext / hPrev / (hPrev + hNext);

      private static double Dc(double hPrev, double hNext)
         => (hNext - hPrev) / hPrev / hNext;

      private static double Dr(double hPrev, double hNext)
         => hPrev / hNext / (hPrev + hNext);
   }
}
