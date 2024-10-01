#include "CLASS.cpp"
#include "Windows.h"
namespace OC {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Сводка для MyForm
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
			double a = 0.0;
			double b = 1.0;
			double epsilon;
			double r;
			unsigned short N;
			unsigned short m = 9;
			double epsilon_granichnoe;
			double eta_0;
			unsigned short global_local_iterations;
			unsigned short max_iterations;
			unsigned short procent_correct;
			unsigned short procent_correct_LNA;
			unsigned short number_of_experiments;
			double correctness;
			Curve_2D = gcnew Hilbert_Curve_2D(a, b, a, b, m, List::Top);
			Curr_start = gcnew CenterSquare(double(), double(), List(), nullptr, nullptr, unsigned short(), double(), double());
			Curr_start = Curve_2D->GetHead();
			Curr_end = gcnew CenterSquare(double(), double(), List(), nullptr, nullptr, unsigned short(), double(), double());
			Curr_end = Curve_2D->GetEnd();
			min_xy = gcnew Extr();
			//
			//TODO: добавьте код конструктора
			//
		}

	protected:
		/// <summary>
		/// Освободить все используемые ресурсы.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::DataVisualization::Charting::Chart^ chart1;
	private: System::Windows::Forms::Button^ button1;



	protected:

	private:
		/// <summary>
		/// Обязательная переменная конструктора.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Требуемый метод для поддержки конструктора — не изменяйте 
		/// содержимое этого метода с помощью редактора кода.
		/// </summary>
		void InitializeComponent(void)
		{
			System::Windows::Forms::DataVisualization::Charting::ChartArea^ chartArea1 = (gcnew System::Windows::Forms::DataVisualization::Charting::ChartArea());
			System::Windows::Forms::DataVisualization::Charting::Legend^ legend1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Legend());
			System::Windows::Forms::DataVisualization::Charting::Series^ series1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Series^ series2 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Title^ title1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Title());
			this->chart1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Chart());
			this->button1 = (gcnew System::Windows::Forms::Button());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->BeginInit();
			this->SuspendLayout();
			// 
			// chart1
			// 
			chartArea1->AxisX->Title = L"Число итераций метода";
			chartArea1->AxisX->TitleFont = (gcnew System::Drawing::Font(L"Tahoma", 14.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			chartArea1->AxisY->Title = L"Доля решенных задач";
			chartArea1->AxisY->TitleFont = (gcnew System::Drawing::Font(L"Tahoma", 14.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			chartArea1->Name = L"ChartArea1";
			this->chart1->ChartAreas->Add(chartArea1);
			legend1->Font = (gcnew System::Drawing::Font(L"Tahoma", 14.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			legend1->IsTextAutoFit = false;
			legend1->Name = L"Legend1";
			this->chart1->Legends->Add(legend1);
			this->chart1->Location = System::Drawing::Point(12, 12);
			this->chart1->Name = L"chart1";
			series1->BorderWidth = 4;
			series1->ChartArea = L"ChartArea1";
			series1->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			series1->Color = System::Drawing::Color::Red;
			series1->Legend = L"Legend1";
			series1->MarkerSize = 10;
			series1->Name = L"AGP";
			series2->BorderWidth = 4;
			series2->ChartArea = L"ChartArea1";
			series2->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			series2->Color = System::Drawing::Color::Blue;
			series2->Legend = L"Legend1";
			series2->MarkerSize = 10;
			series2->Name = L"AGP-LNA";
			this->chart1->Series->Add(series1);
			this->chart1->Series->Add(series2);
			this->chart1->Size = System::Drawing::Size(1642, 739);
			this->chart1->TabIndex = 0;
			this->chart1->Text = L"chart1";
			title1->Font = (gcnew System::Drawing::Font(L"Tahoma", 14.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			title1->Name = L"Title1";
			title1->Text = L"Операционные характеристики АГП";
			this->chart1->Titles->Add(title1);
			// 
			// button1
			// 
			this->button1->Font = (gcnew System::Drawing::Font(L"Tahoma", 14.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button1->Location = System::Drawing::Point(1432, 147);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(158, 36);
			this->button1->TabIndex = 1;
			this->button1->Text = L"ПОСТРОЕНИЕ";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1659, 755);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->chart1);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->EndInit();
			this->ResumeLayout(false);

		}
#pragma endregion
	private:
		double a = 0.0;
		double b = 1.0;
		double epsilon = 0.000005;
		double r = 4.0;
		unsigned short N = 2;
		unsigned short m = 9;
		double epsilon_granichnoe = pow(10.0, -6.0);
		double eta_0 = pow(10.0, -6.0);
		unsigned short global_local_iterations = 10;
		unsigned short max_iterations = 0;
		unsigned short procent_correct;
		unsigned short procent_correct_LNA;
		unsigned short number_of_experiments;
		double correctness = 0.01;
		Hilbert_Curve_2D^ Curve_2D;
		CenterSquare^ Curr_start;
		CenterSquare^ Curr_end;
		Extr^ min_xy;
		System::Void button1_Click(System::Object^ sender, System::EventArgs^ e)
		{
			chart1->ChartAreas[0]->AxisX->Minimum = 0;
			chart1->ChartAreas[0]->AxisX->Maximum = 200;
			chart1->ChartAreas[0]->AxisY->Minimum = 0;
			chart1->ChartAreas[0]->AxisY->Maximum = 100;
			chart1->ChartAreas[0]->AxisX->MajorGrid->Interval = 50;
			chart1->ChartAreas[0]->AxisY->MajorGrid->Interval = 10;
			//
			chart1->Series[0]->Points->AddXY(0, 0);
			chart1->Series[1]->Points->AddXY(0, 0);
			HINSTANCE load_function = LoadLibrary(L"FUNC.dll");
			typedef double (*grsh) (double, double);
			grsh GrishaginFunc = (grsh)GetProcAddress(load_function, "GrishaginFunc");
			while (max_iterations != 200)
			{
				max_iterations += 50;
				procent_correct = 0;
				procent_correct_LNA = 0;
				number_of_experiments = 100;
				while (number_of_experiments != 0)
				{
					Curr_start = Curve_2D->GetHead();
					Curr_end = Curve_2D->GetEnd();
					double minimum = min(GrishaginFunc(Curr_start->Get_x1(), Curr_start->Get_x2()), GrishaginFunc(Curr_end->Get_x1(), Curr_end->Get_x2()));
					Curr_start = Curr_start->GetNext();
					Curr_end = Curr_end->GetPrevious();
					minimum = min((GrishaginFunc(Curr_start->Get_x1(), Curr_start->Get_x2()), GrishaginFunc(Curr_end->Get_x1(), Curr_end->Get_x2())), minimum);
					Curr_start = Curr_start->GetNext();
					Curr_end = Curr_end->GetPrevious();
					while (Curr_start != Curr_end->GetNext())
					{
						minimum = min(GrishaginFunc(Curr_start->Get_x1(), Curr_start->Get_x2()), GrishaginFunc(Curr_end->Get_x1(), Curr_end->Get_x2()));
						Curr_start = Curr_start->GetNext();
						Curr_end = Curr_end->GetPrevious();
						minimum = min((GrishaginFunc(Curr_start->Get_x1(), Curr_start->Get_x2()), GrishaginFunc(Curr_end->Get_x1(), Curr_end->Get_x2())), minimum);
						Curr_start = Curr_start->GetNext();
						Curr_end = Curr_end->GetPrevious();
					}
					min_xy->Base_LNA_1_2_Mer_AGP(a, b, epsilon, r, N, m, epsilon_granichnoe, eta_0, global_local_iterations, max_iterations, false);
					if ((min_xy->GetMin()).Value - minimum <= correctness)
					{
						procent_correct++;
					}
					min_xy->Base_LNA_1_2_Mer_AGP(a, b, epsilon, r, N, m, epsilon_granichnoe, eta_0, global_local_iterations, max_iterations, true);
					if ((min_xy->GetMin()).Value - minimum <= correctness)
					{
						procent_correct_LNA++;
					}
					number_of_experiments--;
				}
				chart1->Series[0]->Points->AddXY(max_iterations, procent_correct);
				chart1->Series[1]->Points->AddXY(max_iterations, procent_correct_LNA);
			}
			FreeLibrary(load_function);
		}
};
}