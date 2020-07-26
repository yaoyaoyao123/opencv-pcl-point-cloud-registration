///*
//���ϵĵ�һ�棬�����鱾��ֱ���Ӳ�ͼ�ؽ���û������ƥ��
//
//���������
//	cam0 = [4152.073 0 1288.147; 0 4152.073 973.571; 0 0 1]
//	cam1 = [4152.073 0 1501.231; 0 4152.073 973.571; 0 0 1]
//	doffs = 213.084
//	baseline = 176.252
//	width = 2872
//	height = 1984
//����ڲ�������
//	K=[fx 0 u0; 0 fy v0; 0 0 1]
//
//	doffs = |u1 - u0|
//*/
//
//#include <pcl/visualization/cloud_viewer.h>
//#include <iostream>  
//#include <pcl/io/io.h>  
//#include <pcl/io/pcd_io.h>  
//#include <opencv2/opencv.hpp>  
//
//using namespace cv;
//using namespace std;
//using namespace pcl;
//
//int user_data;
//// ����ڲ�
//const double u0 = 1288.147;
//const double v0 = 973.571;
//const double fx = 4152.073;
//const double fy = 4152.073;
//const double baseline = 176.252;
//const double doffs = 213.084;	// �����������������x�����ϵĲ��, doffs = |u1 - u0|
//
//void viewerOneOff(visualization::PCLVisualizer& viewer)
//{
//	viewer.setBackgroundColor(0.0, 0.0, 0.0);
//}
//
//int main()
//{
//	// ��������
//	Mat color = imread("dispc0.png"); // RGB
//	Mat depth = imread("disp0.png", IMREAD_UNCHANGED);// depth
//
//	if (color.empty() || depth.empty())
//	{
//		cout << "The image is empty, please check it!" << endl;
//		return -1;
//	}
//
//	// �������ϵ�µĵ���
//	PointCloud<PointXYZRGB>::Ptr cloud(new PointCloud<PointXYZRGB>);
//
//	for (int row = 0; row < depth.rows; row++)
//	{
//		for (int col = 0; col < depth.cols; col++)
//		{
//			ushort d = depth.ptr<ushort>(row)[col];
//
//			if (d == 0)
//				continue;
//			PointXYZRGB p;
//
//			// depth			
// 			p.z = fx * baseline / ((d/16.0) + doffs); // Zc = baseline * f / (d + doffs)
//			p.x = (col - u0) * p.z / fx; // Xc���ң�Yc����Ϊ��
//			p.y = (row - v0) * p.z / fy;
//
//			p.y = -p.y;  // Ϊ������ʾ����x����ά��ת180��
//			p.z = -p.z;
//
//			// RGB
//			p.b = color.ptr<uchar>(row)[col * 3];
//			p.g = color.ptr<uchar>(row)[col * 3 + 1];
//			p.r = color.ptr<uchar>(row)[col * 3 + 2];
//
//			cloud->points.push_back(p);
//		}
//	}
//
//	//cloud->height = depth.rows;
//	//cloud->width = depth.cols;
//	//cloud->points.resize(cloud->height * cloud->width);
//
//	visualization::CloudViewer viewer("Cloud Viewer");
//	viewer.showCloud(cloud);
//	viewer.runOnVisualizationThreadOnce(viewerOneOff);
//
//	while (!viewer.wasStopped())
//	{
//		user_data = 9;
//	}
//	return 0;
//}


///*��֤����SGBM�����õ����Ӳ�ͼ������16�õ���ʵ���Ӳ�ֵ*/
//#include <opencv2/opencv.hpp>
//#include <iostream>
//using namespace std;
//using namespace cv;
//
//void on_mouse(int EVENT, int x, int y, int flags, void* userdata)
//{
//
//	Point p(x, y);
//	switch (EVENT)
//	{
//
//		case EVENT_LBUTTONDOWN:
//		{
//
//		}
//		break;
//
//		case EVENT_RBUTTONDOWN:
//		{
//			cout << "pix coordinate: " << p << endl;
//		}
//
//		break;
//	}
//
//}
//
//int main()
//{
//
//	Mat left = imread("left.jpg", IMREAD_GRAYSCALE);
//	Mat right = imread("right.jpg", IMREAD_GRAYSCALE);
//	Mat disp;
//	int mindisparity = 0;
//	int ndisparities = 16;
//	int SADWindowSize = 11;
//	//SGBM
//	cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(mindisparity, ndisparities, SADWindowSize);
//	int P1 = 8 * left.channels() * SADWindowSize* SADWindowSize;
//	int P2 = 32 * left.channels() * SADWindowSize* SADWindowSize;
//	sgbm->setP1(P1);
//	sgbm->setP2(P2);
//	sgbm->setPreFilterCap(15);
//	sgbm->setUniquenessRatio(10);
//	sgbm->setSpeckleRange(2);
//	sgbm->setSpeckleWindowSize(100);
//	sgbm->setDisp12MaxDiff(1);
//	//sgbm->setMode(cv::StereoSGBM::MODE_HH);
//	sgbm->compute(left, right, disp);
//
//	disp.convertTo(disp, CV_32F, 1.0 / 16);                //����16�õ���ʵ�Ӳ�ֵ
//
//	Mat dispa = imread("dispf.png");
//
//	//cout << disp.ptr<float>(257)[213] << endl;
//	//cout << disp.ptr<float>(42)[310] << endl;
//	//cout << disp.ptr<float>(161)[189] << endl;
//
//	Mat disp8U = Mat(disp.rows, disp.cols, CV_8UC1);       //��ʾ
//	normalize(disp, disp8U, 0, 255, NORM_MINMAX, CV_8UC1);
//
//	//namedWindow("left",WINDOW_AUTOSIZE);
//	//namedWindow("right", WINDOW_AUTOSIZE);
//	//namedWindow("disp", WINDOW_AUTOSIZE);
//	//setMouseCallback("left", on_mouse);
//	//setMouseCallback("right", on_mouse);
//	//setMouseCallback("disp", on_mouse);
//	//imshow("left", left);
//	//imshow("right", right);
//	//imshow("disp", disp8U);
//	waitKey(0);
//	return 0;
//
//}


///*
// �����ƥ����ؽ���άλ��
//*/
//
//#include <iostream>  
//#include <opencv2/opencv.hpp>  
// 
//using namespace cv;
//using namespace std;
//
//// ����ڲ�
//const double u0 = 653.89;
//const double v0 = 349.72;
//const double fx = 706.18;
//const double fy = 706.18;
//const double baseline = 6.17;
//const double doffs = 39.5;	// �����������������x�����ϵĲ��, doffs = |u1 - u0|
//
//// �������
//Point left_point;
//Point right_point;
//
//// ������ά����
//void getClickPointCoord(int event, int x, int y, int flags, void* param);
//void compute3dCoord(Point left_point, Point right_point);
//
//int main()
//{
//
//	Mat left = imread("steroMatchImage/left10.jpg");
//	Mat right = imread("steroMatchImage/right10.jpg");
//
//	/*����У��*/
//	Mat cameraMatrix[2], distCoeffs[2];
//	Mat Q;
//
//	//matlab
//	cameraMatrix[0] = Mat((Mat_<double>(3, 3) << 729.0952, 0.0, 0.0,
//		-0.0867, 729.6895, 0.0,
//		674.2116, 357.1047, 1.0)).t();
//
//	cameraMatrix[1] = Mat((Mat_<double>(3, 3) << 727.3174, 0.0, 0.0,
//		-0.2356, 727.4837, 0.0,
//		634.6699, 342.4619, 1.0)).t();
//
//	//����ϵ��
//	//distCoeffs[0] = (Mat_<double>(4, 1) << 0.0613, 0.1143, -0.0012, 0.0023);
//	//distCoeffs[1] = (Mat_<double>(4, 1) << 0.0806, 0.0180, 0.0014, 0.0011);
//
//	Mat R = Mat((Mat_<double>(3, 3) << 1.0, 1.9618e-04, -0.0022,
//		-1.9726e-04, 1.0, -4.827e-04,
//		0.0022, 4.83e-04, 1.0)).t();
//
//	Mat T = (Mat_<double>(3, 1) << -6.1719, -0.0082, 0.0043); //��������໥������ƽ������
//
//	Mat R1 = (Mat_<double>(3, 3) << 9.9908523344342171e-01, -9.1992112105362279e-03, 4.1762074043410001e-02,
//		9.3234166309280522e-03, 9.9995267073263494e-01, -2.7803231820107406e-03,
//		-4.1734520694857037e-02, 3.1671450510257367e-03, 9.9912371555007995e-01);
//
//	Mat R2 = (Mat_<double>(3, 3) << 9.9913038347395211e-01, -8.7274676713489392e-03, 4.0771413113096469e-02,
//		8.6212419619687286e-03, 9.9995897066935735e-01, 2.7804972418828274e-03,
//		-4.0794006989095910e-02, -2.4265790579512678e-03, 9.9916463143360401e-01);
//
//	Mat Pro1 = (Mat_<double>(3, 4) << 6.1366085042324494e+02, 0., 5.9218479919433594e+02, 0., 0.,
//		6.1366085042324494e+02, 3.3312004470825195e+02, 0., 0., 0., 1.,
//		0.);
//
//	Mat Pro2 = (Mat_<double>(3, 4) << 6.1366085042324494e+02, 0., 5.9218479919433594e+02,
//		-3.7143416039237313e+03, 0., 6.1366085042324494e+02,
//		3.3312004470825195e+02, 0., 0., 0., 1., 0.);
//
//	Rect validRoi[2];
//	Size imageSize = Size(1280, 720);
//
//	//stereoRectify�����ڲκͻ���ϵ�������������������������תR��ƽ�ƾ���T
//	//������ת��ƽ�ƾ���ֽ�Ϊ�����������תһ�����ת����R1��R2�����������������ϵ�µ�ͶӰ����P1��P2����Ȳ���ӳ�����Q
//	//�����õ���bougust����У׼����
//	stereoRectify(cameraMatrix[0], distCoeffs[0],
//		cameraMatrix[1], distCoeffs[1],
//		imageSize, R, T, R1, R2, Pro1, Pro2, Q,
//		CALIB_ZERO_DISPARITY, 1, imageSize, &validRoi[0], &validRoi[1]); //validPixROI1  validPixROI2 - ��ѡ�����������Rect�����ݡ����ڲ����������ض���Ч
//
//	//��ʾ���������ͶӰ����
//	cout << "camL" << cameraMatrix[0] << endl;
//	cout << "camR" << cameraMatrix[1] << endl;
//	cout << "Pro1" << Pro1 <<endl;
//	cout << "Pro2" << Pro2 << endl;
//	cout << "Q" << Q << endl;
//
//	//����������ͼ��У������ӳ���
//	Mat rmap[2][2];
//	initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], R1, Pro1, imageSize, CV_16SC2, rmap[0][0], rmap[0][1]);
//	initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], R2, Pro2, imageSize, CV_16SC2, rmap[1][0], rmap[1][1]);
//
//	//��Դͼ�������ӳ�䣬ӳ�����ΪУ�����ͼ��
//	Mat rleftimg, rrightimg,rcolor;
//	remap(left, rleftimg, rmap[0][0], rmap[0][1], INTER_LINEAR);
//	remap(right, rrightimg, rmap[1][0], rmap[1][1], INTER_LINEAR);
//
//	//�þ��ο�����Ч����
//	Rect vroileft(cvRound(validRoi[0].x), cvRound(validRoi[0].y), cvRound(validRoi[0].width), cvRound(validRoi[0].height));
//	rectangle(rleftimg, vroileft, Scalar(0, 0, 255), 3, 8);
//	Rect vroiright(cvRound(validRoi[1].x), cvRound(validRoi[1].y), cvRound(validRoi[1].width), cvRound(validRoi[1].height));
//	rectangle(rrightimg, vroiright, Scalar(0, 0, 255), 3, 8);
//
//	//�������߽��жԱ�
//	Mat canvas;
//	int w = imageSize.width;
//	int h = imageSize.height;
//	canvas.create(h, w * 2, CV_8UC3);
//	hconcat(rleftimg, rrightimg, canvas);
//
//	for (int j = 0; j < canvas.rows; j += 16)
//		line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
//
//	line(canvas, Point(1280, 0), Point(1280, 720), Scalar(255, 255, 0), 1, 8); // ����ͼ�ָ���
//
//	resize(canvas, canvas, Size(), 0.5, 0.5);
//
//	namedWindow("rectifyResult", WINDOW_AUTOSIZE);
//	setMouseCallback("rectifyResult", getClickPointCoord);
//	imshow("rectifyResult", canvas);
//
//	waitKey(0);
//
//	return 0;
//}
//
//void getClickPointCoord(int event, int x, int y, int flags, void* param)
//{
//	double click_x = x * 2.0;
//	double click_y = y * 2.0;
//
//	if (event == EVENT_LBUTTONDOWN)
//	{
//		// ���е�����ͼ
//		if (click_x <= 1280.0)
//		{
//			left_point = Point(click_x, click_y);
//
//			//cout << "left" << left_point << endl;
//		}
//
//		// ���е�����ͼ
//		if (click_x > 1280.0)
//		{
//			click_x -= 1280.0;
//			right_point = Point(click_x, click_y);
//			//cout << "right" << right_point << endl;
//		}
//
//		// ������ά����
//		compute3dCoord(left_point, right_point);
//	} 
//
//}
//
//void compute3dCoord(Point left_point, Point right_point)
//{
//
//		float d = left_point.x - right_point.x;
//		
//		double z = fx * baseline / d;
//		double x = (left_point.x - u0) * z / fx;
//		double y = (left_point.y - v0) * z / fy;
//
//		cout << "x:" << x << "y:" << y << "z:" << z << endl;
//		//cout << "depth:" << z << endl;
//
//}



///*
//�Լ������������ͼƬ��У����ƥ�䣬�ؽ�����
//*/
//
//#include <pcl/visualization/cloud_viewer.h>
//#include <iostream>  
//#include <pcl/io/io.h>  
//#include <pcl/io/pcd_io.h>  
//#include <opencv2/opencv.hpp>  
//
//using namespace cv;
//using namespace std;
//using namespace pcl;
//
//int user_data;
//// ����ڲ�
//const double u0 = 653.89;
//const double v0 = 349.72;
//const double fx = 706.18;
//const double fy = 706.18;
//const double baseline = 6.17;
//const double doffs = 39.5;	// �����������������x�����ϵĲ��, doffs = |u1 - u0|
//
//void viewerOneOff(visualization::PCLVisualizer& viewer)
//{
//	viewer.setBackgroundColor(0.15, 0.15, 0.15);
//}
//
//int main()
//{
//
//	Mat left = imread("steroMatchImage/left10.jpg", 1);
//	Mat right = imread("steroMatchImage/right10.jpg", 1);
//	Mat color = imread("steroMatchImage/left10.jpg", IMREAD_COLOR);
//
//	Mat canvas_aaa;
//	canvas_aaa.create(left.rows, left.cols * 2, CV_8UC3);
//	hconcat(left, right, canvas_aaa);
//	for (int j = 0; j < canvas_aaa.rows; j += 16)
//		line(canvas_aaa, Point(0, j), Point(canvas_aaa.cols, j), Scalar(0, 255, 0), 1, 8);
//	line(canvas_aaa, Point(1280, 0), Point(1280, 720), Scalar(255, 255, 0), 1, 8);
//	resize(canvas_aaa, canvas_aaa, Size(), 0.5, 0.5);
//	imshow("ssss", canvas_aaa);
//	waitKey(0);
//
//	/*����У��*/
//	Mat cameraMatrix[2], distCoeffs[2];
//	Mat Q;
//
//	//matlab
//	cameraMatrix[0] = Mat((Mat_<double>(3, 3) << 729.0952, 0.0, 0.0,
//		-0.0867, 729.6895, 0.0,
//		674.2116, 357.1047, 1.0)).t();
//
//	cameraMatrix[1] = Mat((Mat_<double>(3, 3) << 727.3174, 0.0, 0.0,
//		-0.2356, 727.4837, 0.0,
//		634.6699, 342.4619, 1.0)).t();
//
//	//����ϵ��
//	//distCoeffs[0] = (Mat_<double>(4, 1) << 0.0613, 0.1143, -0.0012, 0.0023);
//	//distCoeffs[1] = (Mat_<double>(4, 1) << 0.0806, 0.0180, 0.0014, 0.0011);
//
//	Mat R = Mat((Mat_<double>(3, 3) << 1.0, 1.9618e-04, -0.0022,
//		-1.9726e-04, 1.0, -4.827e-04,
//		0.0022, 4.83e-04, 1.0)).t();
//
//	Mat T = (Mat_<double>(3, 1) << -6.1719, -0.0082, 0.0043); //��������໥������ƽ������
//
//	Mat R1 = (Mat_<double>(3, 3) << 9.9908523344342171e-01, -9.1992112105362279e-03, 4.1762074043410001e-02,
//		9.3234166309280522e-03, 9.9995267073263494e-01, -2.7803231820107406e-03,
//		-4.1734520694857037e-02, 3.1671450510257367e-03, 9.9912371555007995e-01);
//
//	Mat R2 = (Mat_<double>(3, 3) << 9.9913038347395211e-01, -8.7274676713489392e-03, 4.0771413113096469e-02,
//		8.6212419619687286e-03, 9.9995897066935735e-01, 2.7804972418828274e-03,
//		-4.0794006989095910e-02, -2.4265790579512678e-03, 9.9916463143360401e-01);
//
//	Mat Pro1 = (Mat_<double>(3, 4) << 6.1366085042324494e+02, 0., 5.9218479919433594e+02, 0., 0.,
//		6.1366085042324494e+02, 3.3312004470825195e+02, 0., 0., 0., 1.,
//		0.);
//
//	Mat Pro2 = (Mat_<double>(3, 4) << 6.1366085042324494e+02, 0., 5.9218479919433594e+02,
//		-3.7143416039237313e+03, 0., 6.1366085042324494e+02,
//		3.3312004470825195e+02, 0., 0., 0., 1., 0.);
//
//	Rect validRoi[2];
//	Size imageSize = Size(1280, 720);
//
//	//stereoRectify�����ڲκͻ���ϵ�������������������������תR��ƽ�ƾ���T
//	//������ת��ƽ�ƾ���ֽ�Ϊ�����������תһ�����ת����R1��R2�����������������ϵ�µ�ͶӰ����P1��P2����Ȳ���ӳ�����Q
//	//�����õ���bougust����У׼����
//	stereoRectify(cameraMatrix[0], distCoeffs[0],
//		cameraMatrix[1], distCoeffs[1],
//		imageSize, R, T, R1, R2, Pro1, Pro2, Q,
//		CALIB_ZERO_DISPARITY, 1, imageSize, &validRoi[0], &validRoi[1]); //validPixROI1  validPixROI2 - ��ѡ�����������Rect�����ݡ����ڲ����������ض���Ч
//
//	//��ʾ���������ͶӰ����
//	cout << "camL" << cameraMatrix[0] << endl;
//	cout << "camR" << cameraMatrix[1] << endl;
//	cout << "Pro1" << Pro1 << endl;
//	cout << "Pro2" << Pro2 << endl;
//	cout << "Q" << Q << endl;
//
//	//����������ͼ��У������ӳ���
//	Mat rmap[2][2];
//	initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], R1, Pro1, imageSize, CV_16SC2, rmap[0][0], rmap[0][1]);
//	initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], R2, Pro2, imageSize, CV_16SC2, rmap[1][0], rmap[1][1]);
//
//	//��Դͼ�������ӳ�䣬ӳ�����ΪУ�����ͼ��
//	Mat rleftimg, rrightimg, rcolor;
//	remap(left, rleftimg, rmap[0][0], rmap[0][1], INTER_LINEAR);
//	remap(right, rrightimg, rmap[1][0], rmap[1][1], INTER_LINEAR);
//	remap(color, rcolor, rmap[0][0], rmap[0][1], INTER_LINEAR);
//
//	Mat sgbmleftimg = rleftimg.clone();
//	Mat sgbmrightimg = rrightimg.clone();
//	//imwrite("rleft.jpg", sgbmleftimg);
//	//imwrite("rright.jpg", sgbmrightimg);
//
//	cvtColor(rleftimg, rleftimg, COLOR_GRAY2BGR);
//	cvtColor(rrightimg, rrightimg, COLOR_GRAY2BGR);
//
//	//�þ��ο�����Ч����
//	Rect vroileft(cvRound(validRoi[0].x), cvRound(validRoi[0].y), cvRound(validRoi[0].width), cvRound(validRoi[0].height));
//	rectangle(rleftimg, vroileft, Scalar(0, 0, 255), 3, 8);
//	Rect vroiright(cvRound(validRoi[1].x), cvRound(validRoi[1].y), cvRound(validRoi[1].width), cvRound(validRoi[1].height));
//	rectangle(rrightimg, vroiright, Scalar(0, 0, 255), 3, 8);
//
//	//�������߽��жԱ�
//	Mat canvas;
//	int w = imageSize.width;
//	int h = imageSize.height;
//	canvas.create(h, w * 2, CV_8UC3);
//	hconcat(rleftimg, rrightimg, canvas);
//
//	for (int j = 0; j < canvas.rows; j += 16)
//		line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
//	line(canvas, Point(1280, 0), Point(1280, 720), Scalar(255, 255, 0), 1, 8);
//
//	resize(canvas, canvas, Size(), 0.5, 0.5);
//	namedWindow("rectifyResult", WINDOW_AUTOSIZE);
//	imshow("rectifyResult", canvas);
//	waitKey(0);
//
//	/*����ƥ��*/
//
//	Mat disp, disp8;
//	int mindisparity = 0;
//	int SADWindowSize = 5;
//
//	//SGBM
//	cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(mindisparity, 256, SADWindowSize);
//
//	int P1 = 8 * sgbmleftimg.channels() * SADWindowSize* SADWindowSize;
//	int P2 = 32 * sgbmleftimg.channels() * SADWindowSize* SADWindowSize;
//
//	sgbm->setP1(P1);
//	sgbm->setP2(P2);
//	sgbm->setBlockSize(5);
//	sgbm->setMinDisparity(0);
//	sgbm->setNumDisparities(256);
//	sgbm->setPreFilterCap(15);
//	sgbm->setUniquenessRatio(10);//��õĴ��۷���ֵ��Ӯ�ˡ��ڶ��õĴ��۷���ֵ�ĸ��ʣ�ͨ������Ϊ5-15֮��Ч���ﵽ��ѣ���Խ����ƥ��Խ��
//	sgbm->setSpeckleRange(1); //�������ص���Ӳ�ֵ������Χ��ͨ������Ϊ1 - 2�ͺ��ˣ����ϵ���ᱻ����16���뵽������
//	sgbm->setSpeckleWindowSize(100);
//	sgbm->setDisp12MaxDiff(1);
//	sgbm->setMode(StereoSGBM::MODE_SGBM);
//	sgbm->compute(sgbmleftimg, sgbmrightimg, disp);
//
//	disp.convertTo(disp, CV_32F, 1.0 / 16.0);                //����16�õ���ʵ�Ӳ�ֵ
//
//	//��ʾ�Ӳ�ͼ
//	Mat disp8U = Mat(disp.rows, disp.cols, CV_8UC1);
//	normalize(disp, disp8U, 0, 255, NORM_MINMAX, CV_8UC1);
//	imshow("disp8U", disp8U);
//	waitKey(0);
//
//	// �������ϵ�µĵ���
//	PointCloud<PointXYZRGB>::Ptr cloud(new PointCloud<PointXYZRGB>);
//	for (int row = 0; row < disp.rows; row++)
//	{
//		for (int col = 0; col < disp.cols; col++)
//		{
//
//			float d = disp.ptr<float>(row)[col];
//			if (d <= 0)
//				continue;
//
//			PointXYZRGB p;
//
//			// depth			
//			p.z = fx * baseline / (d);
//
//			p.x = (col - u0) * p.z / fx;
//			p.y = (row - v0) * p.z / fy;
//
//			p.y = -p.y;  // Ϊ������ʾ����x����ά��ת180��
//			p.z = -p.z;
//
//			// RGB
//			p.b = rcolor.ptr<uchar>(row)[col * 3];
//			p.g = rcolor.ptr<uchar>(row)[col * 3 + 1];
//			p.r = rcolor.ptr<uchar>(row)[col * 3 + 2];
//
//			cloud->points.push_back(p);
//
//		}
//	}
//
//	cloud->height = disp.rows;
//	cloud->width = disp.cols;
//	cloud->points.resize(cloud->height * cloud->width);
//
//	//�������
//	pcl::PCDWriter pcl_writer;
//	pcl_writer.write("pcfile.pcd", *cloud);
//	cout << "save point cloud file success..." << endl;
//
//	//��ʾ����
//	visualization::CloudViewer viewer("Cloud Viewer");
//	viewer.showCloud(cloud);
//	viewer.runOnVisualizationThreadOnce(viewerOneOff);
//
//	while (!viewer.wasStopped())
//	{
//		user_data = 9;
//	}
//	return 0;
//}


///*
//�Լ������������ͼƬ��ֻ��ƥ�䣬�ؽ�����
//*/
//
//#include <pcl/visualization/cloud_viewer.h>
//#include <iostream>  
//#include <pcl/io/io.h>  
//#include <pcl/io/pcd_io.h>  
//#include <opencv2/opencv.hpp>  
//
//using namespace cv;
//using namespace std;
//using namespace pcl;
//
//int user_data;
//
//// ����ڲ�
//const double u0 = 653.89;
//const double v0 = 349.72;
//const double fx = 656.553;
//const double fy = 656.553;
//const double baseline = 59.05;// ����
//const double doffs = 0.0;	// �����������������x�����ϵĲ��, doffs = |u1 - u0|
//
//void viewerOneOff(visualization::PCLVisualizer& viewer)
//{
//	viewer.setBackgroundColor(0.15, 0.15, 0.15);
//}
//
//int main()
//{
//
//	Mat sgbmleftimg = imread("rleft.jpg", 0);
//	Mat sgbmrightimg = imread("rright.jpg", 0);
//	Mat rcolor = imread("rleft.jpg", IMREAD_COLOR);
//
//	/*����ƥ��*/
//
//	Mat disp, disp8;
//	int mindisparity = 0;
//	int SADWindowSize = 5;
//
//	//SGBM
//	cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(mindisparity, 256, SADWindowSize);
//
//	int P1 = 8 * sgbmleftimg.channels() * SADWindowSize* SADWindowSize;
//	int P2 = 32 * sgbmleftimg.channels() * SADWindowSize* SADWindowSize;
//
//	sgbm->setP1(P1);
//	sgbm->setP2(P2);
//	sgbm->setBlockSize(5);
//	sgbm->setMinDisparity(0);
//	sgbm->setNumDisparities(256);
//	sgbm->setPreFilterCap(15);
//	sgbm->setUniquenessRatio(10);//��õĴ��۷���ֵ��Ӯ�ˡ��ڶ��õĴ��۷���ֵ�ĸ��ʣ�ͨ������Ϊ5-15֮��Ч���ﵽ��ѣ���Խ����ƥ��Խ��
//	sgbm->setSpeckleRange(1); //�������ص���Ӳ�ֵ������Χ��ͨ������Ϊ1 - 2�ͺ��ˣ����ϵ���ᱻ����16���뵽������
//	sgbm->setSpeckleWindowSize(100);
//	sgbm->setDisp12MaxDiff(1);
//	sgbm->setMode(StereoSGBM::MODE_SGBM);
//	sgbm->compute(sgbmleftimg, sgbmrightimg, disp);
//
//	disp.convertTo(disp, CV_32F, 1.0 / 16.0);                //����16�õ���ʵ�Ӳ�ֵ
//
//	//��ʾ�Ӳ�ͼ
//	Mat disp8U = Mat(disp.rows, disp.cols, CV_8UC1);
//	normalize(disp, disp8U, 0, 255, NORM_MINMAX, CV_8UC1);
//	imshow("disp8U", disp8U);
//	waitKey(0);
//
//	// �������ϵ�µĵ���
//	PointCloud<PointXYZRGB>::Ptr cloud(new PointCloud<PointXYZRGB>);
//	for (int row = 0; row < disp.rows; row++)
//	{
//		for (int col = 0; col < disp.cols; col++)
//		{
//
//			float d = disp.ptr<float>(row)[col];
//			if (d <= 0)
//				continue;
//
//			PointXYZRGB p;
//
//			// depth			
//			p.z = fx * baseline / (d);
//
//			p.x = (col - u0) * p.z / fx;
//			p.y = (row - v0) * p.z / fy;
//
//			p.y = -p.y;  // Ϊ������ʾ����x����ά��ת180��
//			p.z = -p.z;
//
//			// RGB
//			p.b = rcolor.ptr<uchar>(row)[col * 3];
//			p.g = rcolor.ptr<uchar>(row)[col * 3 + 1];
//			p.r = rcolor.ptr<uchar>(row)[col * 3 + 2];
//
//			cloud->points.push_back(p);
//
//		}
//	}
//
//	cloud->height = disp.rows;
//	cloud->width = disp.cols;
//	cloud->points.resize(cloud->height * cloud->width);
//
//	////�������
//	//pcl::PCDWriter pcl_writer;
//	//pcl_writer.write("pcfile.pcd", *cloud);
//	//cout << "save point cloud file success..." << endl;
//
//	//��ʾ����
//	visualization::CloudViewer viewer("Cloud Viewer");
//	viewer.showCloud(cloud);
//	viewer.runOnVisualizationThreadOnce(viewerOneOff);
//
//	while (!viewer.wasStopped())
//	{
//		user_data = 9;
//	}
//	return 0;
//}


///*
//�Լ������������ͼƬ������Ŀ�������SGBM���ؽ�����
//*/
//
//#include <pcl/visualization/cloud_viewer.h>
//#include <iostream>  
//#include <pcl/io/io.h>  
//#include <pcl/io/pcd_io.h>  
//#include <opencv2/opencv.hpp>  
//
//using namespace cv;
//using namespace std;
//using namespace pcl;
//
//int user_data;
//// ����ڲ�
//const double u0 = 653.89;
//const double v0 = 349.72;
//const double fx = 706.18;
//const double fy = 706.18;
//const double baseline = 6.17;
//const double doffs = 39.5;	// �����������������x�����ϵĲ��, doffs = |u1 - u0|
//
//void viewerOneOff(visualization::PCLVisualizer& viewer)
//{
//	viewer.setBackgroundColor(0.15, 0.15, 0.15);
//}
//
//int main()
//{
//
//	Mat left = imread("steroMatchImage/left10.jpg", 1);
//	Mat right = imread("steroMatchImage/right10.jpg", 1);
//
//	/*����У��*/
//	Mat cameraMatrix[2], distCoeffs[2];
//	Mat Q;
//
//	//matlab
//	cameraMatrix[0] = Mat((Mat_<double>(3, 3) << 729.0952, 0.0, 0.0,
//		-0.0867, 729.6895, 0.0,
//		674.2116, 357.1047, 1.0)).t();
//
//	cameraMatrix[1] = Mat((Mat_<double>(3, 3) << 727.3174, 0.0, 0.0,
//		-0.2356, 727.4837, 0.0,
//		634.6699, 342.4619, 1.0)).t();
//
//	Mat R = Mat((Mat_<double>(3, 3) << 1.0, 1.9618e-04, -0.0022,
//		-1.9726e-04, 1.0, -4.827e-04,
//		0.0022, 4.83e-04, 1.0)).t();
//
//	Mat T = (Mat_<double>(3, 1) << -6.1719, -0.0082, 0.0043); //��������໥������ƽ������
//
//	Mat R1 = (Mat_<double>(3, 3) << 9.9908523344342171e-01, -9.1992112105362279e-03, 4.1762074043410001e-02,
//		9.3234166309280522e-03, 9.9995267073263494e-01, -2.7803231820107406e-03,
//		-4.1734520694857037e-02, 3.1671450510257367e-03, 9.9912371555007995e-01);
//
//	Mat R2 = (Mat_<double>(3, 3) << 9.9913038347395211e-01, -8.7274676713489392e-03, 4.0771413113096469e-02,
//		8.6212419619687286e-03, 9.9995897066935735e-01, 2.7804972418828274e-03,
//		-4.0794006989095910e-02, -2.4265790579512678e-03, 9.9916463143360401e-01);
//
//	Mat Pro1 = (Mat_<double>(3, 4) << 6.1366085042324494e+02, 0., 5.9218479919433594e+02, 0., 0.,
//		6.1366085042324494e+02, 3.3312004470825195e+02, 0., 0., 0., 1.,
//		0.);
//
//	Mat Pro2 = (Mat_<double>(3, 4) << 6.1366085042324494e+02, 0., 5.9218479919433594e+02,
//		-3.7143416039237313e+03, 0., 6.1366085042324494e+02,
//		3.3312004470825195e+02, 0., 0., 0., 1., 0.);
//
//	Rect validRoi[2];
//	Size imageSize = Size(1280, 720);
//
//	//stereoRectify�����ڲκͻ���ϵ�������������������������תR��ƽ�ƾ���T
//	//������ת��ƽ�ƾ���ֽ�Ϊ�����������תһ�����ת����R1��R2�����������������ϵ�µ�ͶӰ����P1��P2����Ȳ���ӳ�����Q
//	//�����õ���bougust����У׼����
//	stereoRectify(cameraMatrix[0], distCoeffs[0],
//		cameraMatrix[1], distCoeffs[1],
//		imageSize, R, T, R1, R2, Pro1, Pro2, Q,
//		CALIB_ZERO_DISPARITY, 1, imageSize, &validRoi[0], &validRoi[1]); //validPixROI1  validPixROI2 - ��ѡ�����������Rect�����ݡ����ڲ����������ض���Ч
//
//	//��ʾ���������ͶӰ����
//	cout << "camL" << cameraMatrix[0] << endl;
//	cout << "camR" << cameraMatrix[1] << endl;
//	cout << "Pro1" << Pro1 << endl;
//	cout << "Pro2" << Pro2 << endl;
//	cout << "Q" << Q << endl;
//
//	//����������ͼ��У������ӳ���
//	Mat rmap[2][2];
//	initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], R1, Pro1, imageSize, CV_16SC2, rmap[0][0], rmap[0][1]);
//	initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], R2, Pro2, imageSize, CV_16SC2, rmap[1][0], rmap[1][1]);
//
//	//��Դͼ�������ӳ�䣬ӳ�����ΪУ�����ͼ��
//	Mat rleftimg, rrightimg;
//	remap(left, rleftimg, rmap[0][0], rmap[0][1], INTER_LINEAR);
//	remap(right, rrightimg, rmap[1][0], rmap[1][1], INTER_LINEAR);
//
//	// ��ȡ����ʾĿ������
//	Mat target_left = rleftimg(Rect(200, 0, 500, 600));
//	Mat target_right = rrightimg(Rect(200, 0, 500, 600));
//	Mat rcolor = target_left;
//
//	imshow("target_left", target_left);
//	imshow("target_right", target_right);
//
//	//�þ��ο�����Ч����
//	Rect vroileft(cvRound(validRoi[0].x), cvRound(validRoi[0].y), cvRound(validRoi[0].width), cvRound(validRoi[0].height));
//	rectangle(rleftimg, vroileft, Scalar(0, 0, 255), 3, 8);
//	Rect vroiright(cvRound(validRoi[1].x), cvRound(validRoi[1].y), cvRound(validRoi[1].width), cvRound(validRoi[1].height));
//	rectangle(rrightimg, vroiright, Scalar(0, 0, 255), 3, 8);
//
//	//�������߽��жԱ�
//	Mat canvas;
//	int w = imageSize.width;
//	int h = imageSize.height;
//	canvas.create(h, w * 2, CV_8UC3);
//	hconcat(rleftimg, rrightimg, canvas);
//
//	for (int j = 0; j < canvas.rows; j += 16)
//		line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
//
//	line(canvas, Point(1280, 0), Point(1280, 720), Scalar(255, 255, 0), 1, 8);
//
//	resize(canvas, canvas, Size(), 0.5, 0.5);
//	//namedWindow("rectifyResult", WINDOW_AUTOSIZE);
//	//imshow("rectifyResult", canvas);
//
//	/*����ƥ��*/
//	Mat disp, disp8;
//	int mindisparity = 0;
//	int SADWindowSize = 5;
//
//	//SGBM
//	cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(mindisparity, 256, SADWindowSize);
//
//	int P1 = 8 * target_left.channels() * SADWindowSize* SADWindowSize;
//	int P2 = 32 * target_left.channels() * SADWindowSize* SADWindowSize;
//
//	sgbm->setP1(P1);
//	sgbm->setP2(P2);
//	sgbm->setBlockSize(5);
//	sgbm->setMinDisparity(0);
//	sgbm->setNumDisparities(256);
//	sgbm->setPreFilterCap(15);
//	sgbm->setUniquenessRatio(10);//��õĴ��۷���ֵ��Ӯ�ˡ��ڶ��õĴ��۷���ֵ�ĸ��ʣ�ͨ������Ϊ5-15֮��Ч���ﵽ��ѣ���Խ����ƥ��Խ��
//	sgbm->setSpeckleRange(1); //�������ص���Ӳ�ֵ������Χ��ͨ������Ϊ1 - 2�ͺ��ˣ����ϵ���ᱻ����16���뵽������
//	sgbm->setSpeckleWindowSize(100);
//	sgbm->setDisp12MaxDiff(1);
//	sgbm->setMode(StereoSGBM::MODE_SGBM);
//	sgbm->compute(target_left, target_right, disp);
//
//	disp.convertTo(disp, CV_32F, 1.0 / 16.0);                //����16�õ���ʵ�Ӳ�ֵ
//
//	//��ʾ�Ӳ�ͼ
//	Mat disp8U = Mat(disp.rows, disp.cols, CV_8UC1);
//	normalize(disp, disp8U, 0, 255, NORM_MINMAX, CV_8UC1);
//	imshow("disp8U", disp8U);
//	waitKey(0);
//
//	// �������ϵ�µĵ���
//	PointCloud<PointXYZRGB>::Ptr cloud(new PointCloud<PointXYZRGB>);
//	for (int row = 0; row < disp.rows; row++)
//	{
//		for (int col = 0; col < disp.cols; col++)
//		{
//
//			float d = disp.ptr<float>(row)[col];
//			if (d <= 0)
//				continue;
//
//			PointXYZRGB p;
//
//			// depth			
//			p.z = fx * baseline / (d);
//
//			p.x = (col - u0) * p.z / fx;
//			p.y = (row - v0) * p.z / fy;
//
//			p.y = -p.y;  // Ϊ������ʾ����x����ά��ת180��
//			p.z = -p.z;
//
//			// RGB
//			p.b = rcolor.ptr<uchar>(row)[col * 3];
//			p.g = rcolor.ptr<uchar>(row)[col * 3 + 1];
//			p.r = rcolor.ptr<uchar>(row)[col * 3 + 2];
//
//			cloud->points.push_back(p);
//
//		}
//	}
//
//	cloud->height = disp.rows;
//	cloud->width = disp.cols;
//	cloud->points.resize(cloud->height * cloud->width);
//
//	visualization::CloudViewer viewer("Cloud Viewer");
//	viewer.showCloud(cloud);
//	viewer.runOnVisualizationThreadOnce(viewerOneOff);
//
//	while (!viewer.wasStopped())
//	{
//		user_data = 9;
//	}
//	return 0;
//}


///*
//���˵��������������ͼƬ������ƥ�䣬�ؽ�����
//
//���������
//	cam0 = [4152.073 0 1288.147; 0 4152.073 973.571; 0 0 1]
//	cam1 = [4152.073 0 1501.231; 0 4152.073 973.571; 0 0 1]
//	doffs = 213.084
//	baseline = 176.252
//	width = 2872
//	height = 1984
//����ڲ�������
//	K=[fx 0 u0; 0 fy v0; 0 0 1]
//
//	doffs = |u1 - u0|
//*/
//
//#include <pcl/visualization/cloud_viewer.h>
//#include <iostream>  
//#include <pcl/io/io.h>  
//#include <pcl/io/pcd_io.h>  
//#include <opencv2/opencv.hpp>  
//
//using namespace cv;
//using namespace std;
//using namespace pcl;
//
//int user_data;
//// ����ڲ�
//const double u0 = 1288.147;
//const double v0 = 973.571;
//const double fx = 4152.073;
//const double fy = 4152.073;
//const double baseline = 176.252;
//const double doffs = 213.084;	// �����������������x�����ϵĲ��, doffs = |u1 - u0|
//
//void viewerOneOff(visualization::PCLVisualizer& viewer)
//{
//	viewer.setBackgroundColor(0.0, 0.0, 0.0);
//}
//
//int main()
//{
//
//	//Mat left = imread("im0.png", 0);
//	//Mat right = imread("im1.png", 0);
//	//Mat color = imread("im0.png", IMREAD_COLOR);
//
//	///*����У��*/
//	//Mat cameraMatrix[2], distCoeffs[2];
//	//Mat Q;
//
//	////matlab
//	//cameraMatrix[0] = (Mat_<double>(3, 3) << 4029.299, 0.0, 1213.198,
//	//	0.0, 4029.299, 975.964,
//	//	0.0, 0.0, 1.0);
//
//	//cameraMatrix[1] = (Mat_<double>(3, 3) << 4029.299, 0.0, 1484.019,
//	//	0.0, 4029.299, 975.964,
//	//	0.0, 0.0, 1.0);
//
//	//Mat R = Mat((Mat_<double>(3, 3) << 1.0, 0.0, 0.0,
//	//	0.0, 1.0, 0.0,
//	//	0.0, 0.0, 1.0)).t();
//
//	//Mat T = (Mat_<double>(3, 1) << -342.789, 0.0, 0.0); //��������໥������ƽ������
//
//	//Mat R1 = (Mat_<double>(3, 3) << 9.9908523344342171e-01, -9.1992112105362279e-03, 4.1762074043410001e-02,
//	//	9.3234166309280522e-03, 9.9995267073263494e-01, -2.7803231820107406e-03,
//	//	-4.1734520694857037e-02, 3.1671450510257367e-03, 9.9912371555007995e-01);
//
//	//Mat R2 = (Mat_<double>(3, 3) << 9.9913038347395211e-01, -8.7274676713489392e-03, 4.0771413113096469e-02,
//	//	8.6212419619687286e-03, 9.9995897066935735e-01, 2.7804972418828274e-03,
//	//	-4.0794006989095910e-02, -2.4265790579512678e-03, 9.9916463143360401e-01);
//
//	//Mat Pro1 = (Mat_<double>(3, 4) << 6.1366085042324494e+02, 0., 5.9218479919433594e+02, 0., 0.,
//	//	6.1366085042324494e+02, 3.3312004470825195e+02, 0., 0., 0., 1.,
//	//	0.);
//
//	//Mat Pro2 = (Mat_<double>(3, 4) << 6.1366085042324494e+02, 0., 5.9218479919433594e+02,
//	//	-3.7143416039237313e+03, 0., 6.1366085042324494e+02,
//	//	3.3312004470825195e+02, 0., 0., 0., 1., 0.);
//
//	//Rect validRoi[2];
//	//Size imageSize = Size(2800, 1908);
//
//	////stereoRectify�����ڲκͻ���ϵ�������������������������תR��ƽ�ƾ���T
//	////������ת��ƽ�ƾ���ֽ�Ϊ�����������תһ�����ת����R1��R2�����������������ϵ�µ�ͶӰ����P1��P2����Ȳ���ӳ�����Q
//	////�����õ���bougust����У׼����
//	//stereoRectify(cameraMatrix[0], distCoeffs[0],
//	//	cameraMatrix[1], distCoeffs[1],
//	//	imageSize, R, T, R1, R2, Pro1, Pro2, Q,
//	//	CALIB_ZERO_DISPARITY, 1, imageSize, &validRoi[0], &validRoi[1]); //validPixROI1  validPixROI2 - ��ѡ�����������Rect�����ݡ����ڲ����������ض���Ч
//
//	////����������ͼ��У������ӳ���
//	//Mat rmap[2][2];
//	//initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], R1, Pro1, imageSize, CV_16SC2, rmap[0][0], rmap[0][1]);
//	//initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], R2, Pro2, imageSize, CV_16SC2, rmap[1][0], rmap[1][1]);
//
//	////��Դͼ�������ӳ�䣬ӳ�����ΪУ�����ͼ��
//	//Mat rleftimg, rrightimg;
//	//remap(left, rleftimg, rmap[0][0], rmap[0][1], INTER_LINEAR);
//	//remap(right, rrightimg, rmap[1][0], rmap[1][1], INTER_LINEAR);
//
//	//Mat sgbmleftimg = rleftimg.clone();
//	//Mat sgbmrightimg = rrightimg.clone();
//
//	//cvtColor(rleftimg, rleftimg, COLOR_GRAY2BGR);
//	//cvtColor(rrightimg, rrightimg, COLOR_GRAY2BGR);
//
//	////�þ��ο�����Ч����
//	//Rect vroileft(cvRound(validRoi[0].x), cvRound(validRoi[0].y), cvRound(validRoi[0].width), cvRound(validRoi[0].height));
//	//rectangle(rleftimg, vroileft, Scalar(0, 0, 255), 3, 8);
//	//Rect vroiright(cvRound(validRoi[1].x), cvRound(validRoi[1].y), cvRound(validRoi[1].width), cvRound(validRoi[1].height));
//	//rectangle(rrightimg, vroiright, Scalar(0, 0, 255), 3, 8);
//
//	////�������߽��жԱ�
//	//Mat canvas;
//	//int w = imageSize.width;
//	//int h = imageSize.height;
//	//canvas.create(h, w * 2, CV_8UC3);
//	//hconcat(rleftimg, rrightimg, canvas);
//
//	//for (int j = 0; j < canvas.rows; j += 16)
//	//	line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
//
//	//line(canvas, Point(2800, 0), Point(2800, 1908), Scalar(255, 255, 0), 1, 8);
//
//	//resize(canvas, canvas, Size(), 0.3, 0.3);
//	//namedWindow("rectifyResult", WINDOW_AUTOSIZE);
//	//imshow("rectifyResult", canvas);
//	//waitKey(0);
//
//	/*����ƥ��*/
//	Mat sgbmleftimg = imread("im0.png", 1);
//	Mat sgbmrightimg = imread("im1.png", 1);
//	Mat color = imread("im0.png", IMREAD_COLOR);
//
//	Mat disp, disp8;
//	int mindisparity = 0;
//	int SADWindowSize = 5;
//
//	//SGBM
//	cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(mindisparity, 240, SADWindowSize);
//
//	int P1 = 8 * sgbmleftimg.channels() * SADWindowSize* SADWindowSize;
//	int P2 = 32 * sgbmleftimg.channels() * SADWindowSize* SADWindowSize;
//
//	sgbm->setP1(P1);
//	sgbm->setP2(P2);
//	sgbm->setBlockSize(5);
//	sgbm->setMinDisparity(0);
//	sgbm->setNumDisparities(240);
//	sgbm->setPreFilterCap(15);
//	sgbm->setUniquenessRatio(10);
//	sgbm->setSpeckleRange(2);
//	sgbm->setSpeckleWindowSize(100);
//	sgbm->setDisp12MaxDiff(1);
//	sgbm->setMode(StereoSGBM::MODE_SGBM);
//
//	sgbm->compute(sgbmleftimg, sgbmrightimg, disp);
//
//	disp.convertTo(disp, CV_32F, 1.0 / 16.0);                //����16�õ���ʵ�Ӳ�ֵ
//	//��ʾ�Ӳ�ͼ
//	Mat disp8U = Mat(disp.rows, disp.cols, CV_8UC1);       
//	normalize(disp, disp8U, 0, 255, NORM_MINMAX, CV_8UC1);
//	resize(disp8U, disp8U, Size(), 0.2, 0.2);
//	imshow("disp8U", disp8U);
//	waitKey(0);
//
//	// �������ϵ�µĵ���
//	PointCloud<PointXYZRGB>::Ptr cloud(new PointCloud<PointXYZRGB>);
//	for (int row = 0; row < disp.rows; row++)
//	{
//		for (int col = 0; col < disp.cols; col++)
//		{
//
//			float d = disp.ptr<float>(row)[col];
//			if (d <= 0)
//				continue;
//
//			PointXYZRGB p;
//
//			// depth			
//			p.z = fx * baseline / (d + doffs); // Zc = baseline * f / (d + doffs)
//			p.x = (col - u0) * p.z / fx; // Xc���ң�Yc����Ϊ��
//			p.y = (row - v0) * p.z / fy;
//
//			p.y = -p.y;  // Ϊ������ʾ����x����ά��ת180��
//			p.z = -p.z;
//
//			// RGB
//			p.b = color.ptr<uchar>(row)[col * 3];
//			p.g = color.ptr<uchar>(row)[col * 3 + 1];
//			p.r = color.ptr<uchar>(row)[col * 3 + 2];
//
//			cloud->points.push_back(p);
//		}
//	}
//
//	visualization::CloudViewer viewer("Cloud Viewer");
//	viewer.showCloud(cloud);
//	viewer.runOnVisualizationThreadOnce(viewerOneOff);
//
//	while (!viewer.wasStopped())
//	{
//		user_data = 9;
//	}
//
//	return 0;
//}


///*
//��������ͼƬ��ֱ���Ӳ�ͼ�ؽ���û������ƥ��
//
//���������
//	cam0 = [4152.073 0 1288.147; 0 4152.073 973.571; 0 0 1]
//	cam1 = [4152.073 0 1501.231; 0 4152.073 973.571; 0 0 1]
//	 doffs = 213.084
//	baseline = 176.252
//	width = 2872
//	height = 1984
//����ڲ�������
//	K=[fx 0 u0; 0 fy v0; 0 0 1]
//
//	doffs = |u1 - u0|
//*/
//
//#include <pcl/visualization/cloud_viewer.h>
//#include <iostream>  
//#include <pcl/io/io.h>  
//#include <pcl/io/pcd_io.h>  
//#include <opencv2/opencv.hpp>  
//#include "PFMReadWrite.h"
//
//using namespace cv;
//using namespace std;
//using namespace pcl;
//
//int user_data;
//// ����ڲ�
//const double u0 = 1213.198;
//const double v0 = 975.964;
//const double fx = 4029.299;
//const double fy = 4029.299;
//const double baseline = 342.789;
//const double doffs = 270.821;	// �����������������x�����ϵĲ��, doffs = |u1 - u0|
//
//void viewerOneOff(visualization::PCLVisualizer& viewer)
//{
//	viewer.setBackgroundColor(0.25, 0.25, 0.25);
//}
//
//int main()
//{
//	// ��������
//	Mat color = imread("dispc0.png"); // RGB
//	Mat depth = imread("disp0.png", IMREAD_UNCHANGED);// depth
//
//	if (color.empty() || depth.empty())
//	{
//		cout << "The image is empty, please check it!" << endl;
//		return -1;
//	}
//
//	// �������ϵ�µĵ���
//	PointCloud<PointXYZRGB>::Ptr cloud(new PointCloud<PointXYZRGB>);
//
//	for (int row = 0; row < depth.rows; row++)
//	{
//		for (int col = 0; col < depth.cols; col++)
//		{
//			ushort d = depth.ptr<ushort>(row)[col];
//
//			if (d == 0)
//				continue;
//			PointXYZRGB p;
//
//			// depth			
// 			p.z = fx * baseline / ((d) + doffs); // Zc = baseline * f / (d + doffs)
//			p.x = (col - u0) * p.z / fx; // Xc���ң�Yc����Ϊ��
//			p.y = (row - v0) * p.z / fy;
//
//			p.y = -p.y;  // Ϊ������ʾ����x����ά��ת180��
//			p.z = -p.z;
//
//			// RGB
//			p.b = color.ptr<uchar>(row)[col * 3];
//			p.g = color.ptr<uchar>(row)[col * 3 + 1];
//			p.r = color.ptr<uchar>(row)[col * 3 + 2];
//
//			cloud->points.push_back(p);
//		}
//	}
//
//	//cloud->height = depth.rows;
//	//cloud->width = depth.cols;
//	//cloud->points.resize(cloud->height * cloud->width);
//
//	visualization::CloudViewer viewer("Cloud Viewer");
//	viewer.showCloud(cloud);
//	viewer.runOnVisualizationThreadOnce(viewerOneOff);
//
//	while (!viewer.wasStopped())
//	{
//		user_data = 9;
//	}
//	return 0;
//}


///* ������׼,��ȡһ���ļ�����תƽ�� */
//#include <pcl/registration/ia_ransac.h>
//#include <pcl/point_types.h>
//#include <pcl/point_cloud.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/features/fpfh.h>
//#include <pcl/search/kdtree.h>
//#include <pcl/io/pcd_io.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/filter.h>
//#include <pcl/registration/icp.h>
//#include <pcl/visualization/pcl_visualizer.h>
//#include <time.h>
//
//using pcl::NormalEstimation;
//using pcl::search::KdTree;
//typedef pcl::PointXYZ PointT;
//typedef pcl::PointCloud<PointT> PointCloud;
//
////���ƿ��ӻ�
//void visualize_pcd(PointCloud::Ptr pcd_src, PointCloud::Ptr pcd_tgt, PointCloud::Ptr pcd_final)
//{
//	//int vp_1, vp_2;
//	// Create a PCLVisualizer object
//	pcl::visualization::PCLVisualizer viewer("registration Viewer");
//	//viewer.createViewPort (0.0, 0, 0.5, 1.0, vp_1);
//   // viewer.createViewPort (0.5, 0, 1.0, 1.0, vp_2);
//
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h(pcd_src, 0, 255, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> tgt_h(pcd_tgt, 255, 0, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> final_h(pcd_final, 0, 0, 255);
//
//	viewer.setBackgroundColor(0.15, 0.15, 0.15, 0); // Setting background to a dark grey
//	viewer.addPointCloud(pcd_src, src_h, "source cloud");
//	viewer.addPointCloud(pcd_tgt, tgt_h, "tgt cloud");
//	viewer.addPointCloud(pcd_final, final_h, "final cloud");
//	//viewer.addCoordinateSystem(1.0);
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(100);
//		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	}
//}
//
////����תƽ�ƾ��������ת�Ƕ�
//void matrix2angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle)
//{
//	double ax, ay, az;
//	if (result_trans(2, 0) == 1 || result_trans(2, 0) == -1)
//	{
//		az = 0;
//		double dlta;
//		dlta = atan2(result_trans(0, 1), result_trans(0, 2));
//		if (result_trans(2, 0) == -1)
//		{
//			ay = M_PI / 2;
//			ax = az + dlta;
//		}
//		else
//		{
//			ay = -M_PI / 2;
//			ax = -az + dlta;
//		}
//	}
//	else
//	{
//		ay = -asin(result_trans(2, 0));
//		ax = atan2(result_trans(2, 1) / cos(ay), result_trans(2, 2) / cos(ay));
//		az = atan2(result_trans(1, 0) / cos(ay), result_trans(0, 0) / cos(ay));
//	}
//	result_angle << ax, ay, az;
//}
//
////����תƽ�ƾ������ƽ�ƾ���
//void matrix2translate(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_translate)
//{
//	result_translate << result_trans(0, 3), result_trans(1, 3), result_trans(2, 3);
//} 
//
//
//
//int main(int argc, char** argv)
//{
//	// ����Դ�����ļ�
//	PointCloud::Ptr cloud_src_o(new PointCloud);
//	pcl::io::loadPCDFile("rabbit_gra.pcd", *cloud_src_o);
//
//	clock_t start = clock();
//
//	// ȥ��Դ����NAN��
//	std::vector<int> indices_src;
//	pcl::removeNaNFromPointCloud(*cloud_src_o, *cloud_src_o, indices_src);
//	std::cout << "remove *cloud_src nan" << endl;
//	// Դ�����²���
//	pcl::VoxelGrid<pcl::PointXYZ> voxel_grid_2;
//	float leaf_size_src = 1.0;
//	voxel_grid_2.setLeafSize(leaf_size_src, leaf_size_src, leaf_size_src);
//	voxel_grid_2.setInputCloud(cloud_src_o);
//	PointCloud::Ptr cloud_src(new PointCloud);// �²�����ĵ���
//	voxel_grid_2.filter(*cloud_src);
//	std::cout << "down size *cloud_src_o.pcd from " << cloud_src_o->size() << "to" << cloud_src->size() << endl;
//
//	//��תƽ������Ŀ�����
//	Eigen::Affine3f transform = Eigen::Affine3f::Identity();
//	Eigen::Vector3f ANGLE_origin;
//	Eigen::Vector3f TRANSLATE_origin;
//	ANGLE_origin << M_PI / 6, M_PI / 8, -M_PI / 3;
//	TRANSLATE_origin << 20.0, -10.0, 6.0;
//
//	transform.translation() << TRANSLATE_origin;
//	transform.rotate(Eigen::AngleAxisf(ANGLE_origin.z(), Eigen::Vector3f::UnitZ()));
//	transform.rotate(Eigen::AngleAxisf(ANGLE_origin.x(), Eigen::Vector3f::UnitX()));
//	transform.rotate(Eigen::AngleAxisf(ANGLE_origin.y(), Eigen::Vector3f::UnitY()));
//	std::cout << "transform matrix " << "\n" << transform.matrix() << endl;
//
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tar(new pcl::PointCloud<pcl::PointXYZ>());
//	pcl::transformPointCloud(*cloud_src, *cloud_tar, transform);
//
//	// ������淨��
//	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne_src; // Դ
//	ne_src.setInputCloud(cloud_src);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_src.setSearchMethod(tree_src);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_src.setRadiusSearch(0.02);
//	ne_src.compute(*cloud_src_normals);
//
//	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne_tar; // Ŀ��
//	ne_tar.setInputCloud(cloud_tar);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_tar(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_tar.setSearchMethod(tree_tar);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_tar_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_tar.setRadiusSearch(0.02);
//	ne_tar.compute(*cloud_tar_normals);
//
//	//����FPFH
//	pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_src;
//	fpfh_src.setInputCloud(cloud_src);
//	fpfh_src.setInputNormals(cloud_src_normals);
//	pcl::search::KdTree<PointT>::Ptr tree_src_fpfh(new pcl::search::KdTree<PointT>);
//	fpfh_src.setSearchMethod(tree_src_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_src(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_src.setRadiusSearch(0.05);
//	fpfh_src.compute(*fpfhs_src);
//	std::cout << "compute *cloud_src fpfh" << endl;
//
//	pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_tar;
//	fpfh_tar.setInputCloud(cloud_tar);
//	fpfh_tar.setInputNormals(cloud_tar_normals);
//	pcl::search::KdTree<PointT>::Ptr tree_tgt_fpfh(new pcl::search::KdTree<PointT>);
//	fpfh_tar.setSearchMethod(tree_tgt_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_tar(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_tar.setRadiusSearch(0.05);
//	fpfh_tar.compute(*fpfhs_tar);
//	std::cout << "compute *cloud_tar fpfh" << endl;
//
//	//SAC��׼
//	pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> scia;
//	scia.setInputSource(cloud_src);
//	scia.setInputTarget(cloud_tar);
//	scia.setSourceFeatures(fpfhs_src);
//	scia.setTargetFeatures(fpfhs_tar);
//	//scia.setMinSampleDistance(1);
//	//scia.setNumberOfSamples(2);
//	//scia.setCorrespondenceRandomness(20);
//	PointCloud::Ptr sac_result(new PointCloud);
//	scia.align(*sac_result);
//	std::cout << "sac has converged:" << scia.hasConverged() << "  score: " << scia.getFitnessScore() << endl;
//	Eigen::Matrix4f sac_trans;
//	sac_trans = scia.getFinalTransformation();
//	std::cout << sac_trans << endl;
//	clock_t sac_time = clock();
//
//	//icp��׼
//	PointCloud::Ptr icp_result(new PointCloud);
//	pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
//	icp.setInputSource(cloud_src);
//	icp.setInputTarget(cloud_tar);
//	// ��������Ӧ����� (e.g., correspondences with higher distances will be ignored)
//	icp.setMaxCorrespondenceDistance(40);
//	// ����������
//	icp.setMaximumIterations(100);
//	// ���α仯����֮��Ĳ�ֵ
//	icp.setTransformationEpsilon(1e-10);
//	// �������
//	icp.setEuclideanFitnessEpsilon(0.2);
//	icp.align(*icp_result, sac_trans);
//
//	clock_t end = clock();
//	cout << "total time: " << (double)(end - start) / (double)CLOCKS_PER_SEC << " s" << endl;
//	//���㷨�ߺ͵�����ֱ��ͼ��ʱ��Ҳ����SAC����
//	cout << "sac time: " << (double)(sac_time - start) / (double)CLOCKS_PER_SEC << " s" << endl;
//	cout << "icp time: " << (double)(end - sac_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	std::cout << "ICP has converged:" << icp.hasConverged()
//		<< " score: " << icp.getFitnessScore() << std::endl;
//	Eigen::Matrix4f icp_trans;
//	icp_trans = icp.getFinalTransformation();
//
//	//cout<<"ransformationProbability"<<icp.getTransformationProbability()<<endl;
//	std::cout << icp_trans << endl;
//	//ʹ�ü���ı任���б任
//	pcl::transformPointCloud(*cloud_src, *icp_result, icp_trans);
//
//	//�������
//	Eigen::Vector3f ANGLE_result;
//	Eigen::Vector3f TRANSLATE_result;
//	matrix2angle(icp_trans, ANGLE_result);
//	matrix2translate(icp_trans, TRANSLATE_result);
//
//	double A_error_x = fabs(ANGLE_result(0)) - fabs(ANGLE_origin(0));
//	double A_error_y = fabs(ANGLE_result(1)) - fabs(ANGLE_origin(1));
//	double A_error_z = fabs(ANGLE_result(2)) - fabs(ANGLE_origin(2));
//	double T_error_x = fabs(TRANSLATE_result(0)) - fabs(TRANSLATE_origin(0));
//	double T_error_y = fabs(TRANSLATE_result(1)) - fabs(TRANSLATE_origin(1));
//	double T_error_z = fabs(TRANSLATE_result(2)) - fabs(TRANSLATE_origin(2));
//
//	cout << "original angle in x y z:\n" << ANGLE_origin << endl;
//	cout << "result angle in x y z:\n" << ANGLE_result << endl;
//	cout << "angle error in aixs_x: " << A_error_x << "  angle error in aixs_y: " << A_error_y << "  angle error in aixs_z: " << A_error_z << endl;
//	cout << "original translate in x y z:\n" << TRANSLATE_origin << endl;
//	cout << "result translate in x y z:\n" << TRANSLATE_result << endl;
//	cout << "translate error in aixs_x: " << T_error_x << "  translate error in aixs_y: " << T_error_y << "  translate error in aixs_z: " << T_error_z << endl;
//
//	//���ӻ�����ɫ��Դ���ƣ���ɫ��Ŀ����ƣ���ɫ����׼֮��ĵ���
//	visualize_pcd(cloud_src, cloud_tar, icp_result);
//	return (0);
//}


///* ��ͬ�ܶȵĵ��ƽ�����׼ */
//#include <pcl/registration/ia_ransac.h>
//#include <pcl/point_types.h>
//#include <pcl/point_cloud.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/features/fpfh.h>
//#include <pcl/search/kdtree.h>
//#include <pcl/io/pcd_io.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/filter.h>
//#include <pcl/registration/icp.h>
//#include <pcl/visualization/pcl_visualizer.h>
//#include <time.h>
//
//using pcl::NormalEstimation;
//using pcl::search::KdTree;
//typedef pcl::PointXYZ PointT;
//typedef pcl::PointCloud<PointT> PointCloud;
//
////���ƿ��ӻ�
//void visualize_pcd(PointCloud::Ptr pcd_src, PointCloud::Ptr pcd_tgt, PointCloud::Ptr pcd_final)
//{
//	pcl::visualization::PCLVisualizer viewer("registration Viewer");
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h(pcd_src, 0, 255, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> tgt_h(pcd_tgt, 255, 0, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> final_h(pcd_final, 0, 0, 255);
//
//	viewer.setBackgroundColor(0.15, 0.15, 0.15, 0); // Setting background to a dark grey
//	viewer.addPointCloud(pcd_src, src_h, "source cloud");
//	viewer.addPointCloud(pcd_tgt, tgt_h, "tgt cloud");
//	viewer.addPointCloud(pcd_final, final_h, "final cloud");
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(100);
//		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	}
//}
//
////����תƽ�ƾ��������ת�Ƕ�
//void matrix2angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle)
//{
//	double ax, ay, az;
//	if (result_trans(2, 0) == 1 || result_trans(2, 0) == -1)
//	{
//		az = 0;
//		double dlta;
//		dlta = atan2(result_trans(0, 1), result_trans(0, 2));
//		if (result_trans(2, 0) == -1)
//		{
//			ay = M_PI / 2;
//			ax = az + dlta;
//		}
//		else
//		{
//			ay = -M_PI / 2;
//			ax = -az + dlta;
//		}
//	}
//	else
//	{
//		ay = -asin(result_trans(2, 0));
//		ax = atan2(result_trans(2, 1) / cos(ay), result_trans(2, 2) / cos(ay));
//		az = atan2(result_trans(1, 0) / cos(ay), result_trans(0, 0) / cos(ay));
//	}
//	result_angle << ax, ay, az;
//}
//
////����תƽ�ƾ������ƽ�ƾ���
//void matrix2translate(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_translate)
//{
//	result_translate << result_trans(0, 3), result_trans(1, 3), result_trans(2, 3);
//}
//
//
//
//int main(int argc, char** argv)
//{
//	// ����Դ�����ļ�
//	PointCloud::Ptr cloud_src_o(new PointCloud);
//	pcl::io::loadPCDFile("rabbit_gra.pcd", *cloud_src_o);
//
//	clock_t start = clock();
//
//	// ȥ��Դ����NAN��
//	std::vector<int> indices_src;
//	pcl::removeNaNFromPointCloud(*cloud_src_o, *cloud_src_o, indices_src);
//	std::cout << "remove *cloud_src nan" << endl;
//
//	// Դ�����²���
//	pcl::VoxelGrid<pcl::PointXYZ> voxel_grid_2;
//	float leaf_size_src = 1.0;
//	voxel_grid_2.setLeafSize(leaf_size_src, leaf_size_src, leaf_size_src);
//	voxel_grid_2.setInputCloud(cloud_src_o);
//	PointCloud::Ptr cloud_src_dowm(new PointCloud);// �²�����ĵ���
//	voxel_grid_2.filter(*cloud_src_dowm);
//	std::cout << "down size *cloud_src_o.pcd from " << cloud_src_o->size() << "to" << cloud_src_dowm->size() << endl;
//
//	// ��תƽ������Ŀ�����
//	Eigen::Affine3f transform = Eigen::Affine3f::Identity();
//	Eigen::Vector3f ANGLE_origin;
//	Eigen::Vector3f TRANSLATE_origin;
//	ANGLE_origin << M_PI / 6, M_PI / 8, -M_PI / 3;
//	TRANSLATE_origin << 20.0, -10.0, 6.0;
//
//	transform.translation() << TRANSLATE_origin;
//	transform.rotate(Eigen::AngleAxisf(ANGLE_origin.z(), Eigen::Vector3f::UnitZ()));
//	transform.rotate(Eigen::AngleAxisf(ANGLE_origin.x(), Eigen::Vector3f::UnitX()));
//	transform.rotate(Eigen::AngleAxisf(ANGLE_origin.y(), Eigen::Vector3f::UnitY()));
//	std::cout << "transform matrix " << "\n" << transform.matrix() << endl;
//
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tar(new pcl::PointCloud<pcl::PointXYZ>());
//	pcl::transformPointCloud(*cloud_src_dowm, *cloud_tar, transform);
//
//	// ������淨��
//	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne_src; // Դ
//	ne_src.setInputCloud(cloud_src_o);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_src.setSearchMethod(tree_src);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_src.setRadiusSearch(0.02);
//	ne_src.compute(*cloud_src_normals);
//
//	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne_tar; // Ŀ��
//	ne_tar.setInputCloud(cloud_tar);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_tar(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_tar.setSearchMethod(tree_tar);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_tar_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_tar.setRadiusSearch(0.02);
//	ne_tar.compute(*cloud_tar_normals);
//
//	//����FPFH
//	pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_src;
//	fpfh_src.setInputCloud(cloud_src_o);
//	fpfh_src.setInputNormals(cloud_src_normals);
//	pcl::search::KdTree<PointT>::Ptr tree_src_fpfh(new pcl::search::KdTree<PointT>);
//	fpfh_src.setSearchMethod(tree_src_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_src(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_src.setRadiusSearch(0.05);
//	fpfh_src.compute(*fpfhs_src);
//	std::cout << "compute *cloud_src fpfh" << endl;
//
//	pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_tar;
//	fpfh_tar.setInputCloud(cloud_tar);
//	fpfh_tar.setInputNormals(cloud_tar_normals);
//	pcl::search::KdTree<PointT>::Ptr tree_tgt_fpfh(new pcl::search::KdTree<PointT>);
//	fpfh_tar.setSearchMethod(tree_tgt_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_tar(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_tar.setRadiusSearch(0.05);
//	fpfh_tar.compute(*fpfhs_tar);
//	std::cout << "compute *cloud_tar fpfh" << endl;
//
//	//SAC��׼
//	pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> scia;
//	scia.setInputSource(cloud_src_o);
//	scia.setInputTarget(cloud_tar);
//	scia.setSourceFeatures(fpfhs_src);
//	scia.setTargetFeatures(fpfhs_tar);
//	//scia.setMinSampleDistance(1);
//	//scia.setNumberOfSamples(2);
//	//scia.setCorrespondenceRandomness(20);
//	PointCloud::Ptr sac_result(new PointCloud);
//	scia.align(*sac_result);
//	std::cout << "sac has converged:" << scia.hasConverged() << "  score: " << scia.getFitnessScore() << endl;
//	Eigen::Matrix4f sac_trans;
//	sac_trans = scia.getFinalTransformation();
//	std::cout << sac_trans << endl;
//	clock_t sac_time = clock();
//
//	//icp��׼
//	PointCloud::Ptr icp_result(new PointCloud);
//	pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
//	icp.setInputSource(cloud_src_o);
//	icp.setInputTarget(cloud_tar);
//	// ��������Ӧ����� (e.g., correspondences with higher distances will be ignored)
//	icp.setMaxCorrespondenceDistance(40);
//	// ����������
//	icp.setMaximumIterations(100);
//	// ���α仯����֮��Ĳ�ֵ
//	icp.setTransformationEpsilon(1e-10);
//	// �������
//	icp.setEuclideanFitnessEpsilon(0.2);
//	icp.align(*icp_result, sac_trans);
//
//	clock_t end = clock();
//	cout << "total time: " << (double)(end - start) / (double)CLOCKS_PER_SEC << " s" << endl;
//	//���㷨�ߺ͵�����ֱ��ͼ��ʱ��Ҳ����SAC����
//	cout << "sac time: " << (double)(sac_time - start) / (double)CLOCKS_PER_SEC << " s" << endl;
//	cout << "icp time: " << (double)(end - sac_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	std::cout << "ICP has converged:" << icp.hasConverged()
//		<< " score: " << icp.getFitnessScore() << std::endl;
//	Eigen::Matrix4f icp_trans;
//	icp_trans = icp.getFinalTransformation();
//
//	//cout<<"ransformationProbability"<<icp.getTransformationProbability()<<endl;
//	std::cout << icp_trans << endl;
//	//ʹ�ü���ı任���б任
//	pcl::transformPointCloud(*cloud_src_o, *icp_result, icp_trans);
//
//	//�������
//	Eigen::Vector3f ANGLE_result;
//	Eigen::Vector3f TRANSLATE_result;
//	matrix2angle(icp_trans, ANGLE_result);
//	matrix2translate(icp_trans, TRANSLATE_result);
//
//	double A_error_x = fabs(ANGLE_result(0)) - fabs(ANGLE_origin(0));
//	double A_error_y = fabs(ANGLE_result(1)) - fabs(ANGLE_origin(1));
//	double A_error_z = fabs(ANGLE_result(2)) - fabs(ANGLE_origin(2));
//	double T_error_x = fabs(TRANSLATE_result(0)) - fabs(TRANSLATE_origin(0));
//	double T_error_y = fabs(TRANSLATE_result(1)) - fabs(TRANSLATE_origin(1));
//	double T_error_z = fabs(TRANSLATE_result(2)) - fabs(TRANSLATE_origin(2));
//
//	cout << "original angle in x y z:\n" << ANGLE_origin << endl;
//	cout << "result angle in x y z:\n" << ANGLE_result << endl;
//	cout << "angle error in aixs_x: " << A_error_x << "  angle error in aixs_y: " << A_error_y << "  angle error in aixs_z: " << A_error_z << endl;
//	cout << "original translate in x y z:\n" << TRANSLATE_origin << endl;
//	cout << "result translate in x y z:\n" << TRANSLATE_result << endl;
//	cout << "translate error in aixs_x: " << T_error_x << "  translate error in aixs_y: " << T_error_y << "  translate error in aixs_z: " << T_error_z << endl;
//
//	//���ӻ�����ɫ��Դ���ƣ���ɫ��Ŀ����ƣ���ɫ����׼֮��ĵ���
//	visualize_pcd(cloud_src_o, cloud_tar, icp_result);
//	return (0);
//}


///* ��ȡT����ĸ��stl�ļ������е�����׼ */
//#include <iostream>
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/common/io.h>
//#include <pcl/point_cloud.h>
//
//#include <pcl/PolygonMesh.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <pcl/io/vtk_lib_io.h>
//#include <pcl/io/io.h>
//#include <pcl/io/pcd_io.h>//pcd ��д����ص�ͷ�ļ���
//#include <pcl/point_types.h> //PCL��֧�ֵĵ�����ͷ�ļ�
//#include <pcl/registration/ia_ransac.h>
//#include <pcl/point_cloud.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/features/normal_3d_omp.h>
//#include <pcl/features/fpfh.h>
//#include <pcl/features/fpfh_omp.h>
//#include <pcl/search/kdtree.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/filter.h>
//#include <pcl/registration/icp.h>
//#include <pcl/visualization/pcl_visualizer.h>
//#include <time.h>
//
//int user_data;
//using namespace std;
//
//void viewerOneOff(pcl::visualization::PCLVisualizer& viewer) {
//	viewer.setBackgroundColor(0.15, 0.15, 0.15);   //���ñ�����ɫ
//}
//
//// �߿����
////pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame(new pcl::PointCloud<pcl::PointXYZ>());
//
//void createTluomuWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame);
//void createPartTluomuWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame);
//
//int main()
//{
//
///* ��ʾmesh */
//
//	//// ������ʾ����
//	//pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
//	//pcl::PolygonMesh mesh;
//	//pcl::io::loadPolygonFileSTL("./jixielingjian/Tluomu.stl", mesh);
//	//pcl::io::loadPCDFile<pcl::PointXYZ>("./jixielingjian/Tluomu.pcd", *cloud);
//	//pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame(new pcl::PointCloud<pcl::PointXYZ>());
//	//pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_part_wire_frame(new pcl::PointCloud<pcl::PointXYZ>());
//	//createTluomuWireFramePointCloud(cloud_wire_frame);
//	//createPartTluomuWireFramePointCloud(cloud_part_wire_frame);
//
//	//// �����ӿ�
//	//boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("The name"));						  //�������ӻ�viewer����
//	//int v1(0), v2(1);
//	//viewer->createViewPort(0.0, 0.0, 0.5, 1, v1);//xmin,ymin,xmax,ymax
//	//viewer->createViewPort(0.5, 0, 1, 1, v2);
//	//viewer->setBackgroundColor(0.15, 0.15, 0.15, v1);
//	//viewer->setBackgroundColor(0.15, 0.15, 0.15, v2); 
//
//	////viewer->addPolygonMesh(mesh, "mesh", v1); // ��ʾmesh
//	//viewer->addPointCloud(cloud_wire_frame, "wf_points", v1); // ��ʾģ���߿�cloud
//	//viewer->addPointCloud(cloud_part_wire_frame, "points", v2); // ��ʾģ��cloud
//	////viewer->setRepresentationToSurfaceForAllActors(); //����ģ������Ƭ��ʽ��ʾ  
//	////viewer->setRepresentationToPointsForAllActors(); //����ģ���Ե���ʽ��ʾ  
//	////viewer->setRepresentationToWireframeForAllActors();  //����ģ�����߿�ͼģʽ��ʾ
//	//viewer->addCoordinateSystem(20.0); // �������ָʾ�ᣬ���ű���
//	//while (!viewer->wasStopped())
//	//{
//	//	viewer->spinOnce(100);
//	//	boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	//}
//	//return 0;
//
//
//
/////* ������ƺ��߿������׼ */
//	//// ������ʾ����
//	//pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
//	//pcl::io::loadPCDFile<pcl::PointXYZ>("./jixielingjian/Tluomu.pcd", *cloud);
//	////pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame(new pcl::PointCloud<pcl::PointXYZ>());
//	//pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_part_wire_frame(new pcl::PointCloud<pcl::PointXYZ>());
//	////createTluomuWireFramePointCloud(cloud_wire_frame);
//	//createPartTluomuWireFramePointCloud(cloud_part_wire_frame);
//
//	//// ��תƽ������Դ����
//	//Eigen::Affine3f transform = Eigen::Affine3f::Identity();
//	//Eigen::Vector3f ANGLE_origin;
//	//Eigen::Vector3f TRANSLATE_origin;
//	//ANGLE_origin << M_PI / 6, M_PI / 8, -M_PI / 3;
//	//TRANSLATE_origin << 50.0, -30.0, 70.0;
//	//transform.translation() << TRANSLATE_origin;
//	//transform.rotate(Eigen::AngleAxisf(ANGLE_origin.z(), Eigen::Vector3f::UnitZ()));
//	//transform.rotate(Eigen::AngleAxisf(ANGLE_origin.x(), Eigen::Vector3f::UnitX()));
//	//transform.rotate(Eigen::AngleAxisf(ANGLE_origin.y(), Eigen::Vector3f::UnitY()));
//	//pcl::transformPointCloud(*cloud_part_wire_frame, *cloud_part_wire_frame, transform);
//
//	//// ������淨��
//	//clock_t fpfh_start_time = clock();
//	//pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_src; // Դ
//	//ne_src.setInputCloud(cloud_part_wire_frame);
//	//pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree< pcl::PointXYZ>());
//	//ne_src.setSearchMethod(tree_src);
//	//pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals(new pcl::PointCloud< pcl::Normal>);
//	//ne_src.setRadiusSearch(2.0);
//	//ne_src.setNumberOfThreads(8);
//	//ne_src.compute(*cloud_src_normals);
//
//	//pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_tar; // Ŀ��
//	//ne_tar.setInputCloud(cloud);
//	//pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_tar(new pcl::search::KdTree< pcl::PointXYZ>());
//	//ne_tar.setSearchMethod(tree_tar);
//	//pcl::PointCloud<pcl::Normal>::Ptr cloud_tar_normals(new pcl::PointCloud< pcl::Normal>);
//	//ne_tar.setRadiusSearch(2.0);
//	//ne_src.setNumberOfThreads(8);
//	//ne_tar.compute(*cloud_tar_normals);
//
//	////����FPFH
//	//pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_src;
//	//fpfh_src.setInputCloud(cloud_part_wire_frame);
//	//fpfh_src.setInputNormals(cloud_src_normals);
//	//pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_src_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	//fpfh_src.setSearchMethod(tree_src_fpfh);
//	//pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_src(new pcl::PointCloud<pcl::FPFHSignature33>());
//	//fpfh_src.setRadiusSearch(2.0);
//	//fpfh_src.setNumberOfThreads(8);
//	//fpfh_src.compute(*fpfhs_src);
//	//std::cout << "compute *cloud_src fpfh" << endl;
//
//	//pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_tar;
//	//fpfh_tar.setInputCloud(cloud);
//	//fpfh_tar.setInputNormals(cloud_tar_normals);
//	//pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_tgt_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	//fpfh_tar.setSearchMethod(tree_tgt_fpfh);
//	//pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_tar(new pcl::PointCloud<pcl::FPFHSignature33>());
//	//fpfh_tar.setRadiusSearch(0.05);
//	//fpfh_src.setNumberOfThreads(8);
//	//fpfh_tar.compute(*fpfhs_tar);
//
//	//std::cout << "compute *cloud_tar fpfh" << endl;
//	//clock_t fpfh_end_time = clock();
//	//cout << "fpfh time: " << (double)(fpfh_end_time - fpfh_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	////SAC��׼
//	//clock_t sac_start_time = clock();
//
//	//pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> scia;
//	//scia.setInputSource(cloud_part_wire_frame);
//	//scia.setInputTarget(cloud);
//	//scia.setSourceFeatures(fpfhs_src);
//	//scia.setTargetFeatures(fpfhs_tar);
//	//pcl::PointCloud<pcl::PointXYZ>::Ptr sac_result(new pcl::PointCloud<pcl::PointXYZ>);
//	//scia.align(*sac_result);
//	//Eigen::Matrix4f sac_trans = scia.getFinalTransformation();
//	//clock_t sac_end_time = clock();
//	//cout << "sac time: " << (double)(sac_end_time - sac_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	////icp��׼
//	//clock_t icp_start_time = clock();
//	//pcl::PointCloud<pcl::PointXYZ>::Ptr icp_result(new pcl::PointCloud<pcl::PointXYZ>);
//	//pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
//	//icp.setInputSource(cloud_part_wire_frame);
//	//icp.setInputTarget(cloud);
//	//// ��������Ӧ����� (e.g., correspondences with higher distances will be ignored)
//	//icp.setMaxCorrespondenceDistance(40);
//	//// ����������
//	//icp.setMaximumIterations(100);
//	//// ���α仯����֮��Ĳ�ֵ
//	//icp.setTransformationEpsilon(1e-10);
//	//// �������
//	//icp.setEuclideanFitnessEpsilon(0.2);
//	//icp.align(*icp_result, sac_trans);
//	//clock_t icp_end_time = clock();
//	//cout << "icp time: " << (double)(icp_end_time - icp_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//// ��ʾ����
//	//// ��ɫ��Դ����
//	//// ��ɫ��Ŀ�����
//	//// ��ɫ����׼����
//	//pcl::visualization::PCLVisualizer viewer("registration Viewer");
//	//pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h(cloud_part_wire_frame, 0, 255, 0);
//	//pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> tgt_h(cloud, 255, 0, 0);
//	//pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> final_h(icp_result, 0, 0, 255);
//	//viewer.setBackgroundColor(0.15, 0.15, 0.15);
//
//	//viewer.addPointCloud(cloud_part_wire_frame, src_h, "tar_points");
//	//viewer.addPointCloud(cloud, tgt_h, "src_points");
//	//viewer.addPointCloud(icp_result, final_h, "icp_points");
//
//	//viewer.addCoordinateSystem(20.0); // �������ָʾ�ᣬ���ű���
//	//while (!viewer.wasStopped())
//	//{
//	//	viewer.spinOnce(100);
//	//	boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	//}
//	//return 0;
//
//
//
//
///* �߿���߿������׼ */
//
//	// ������ʾ����
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame(new pcl::PointCloud<pcl::PointXYZ>());
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_part_wire_frame(new pcl::PointCloud<pcl::PointXYZ>());
//	createTluomuWireFramePointCloud(cloud_wire_frame);
//	createPartTluomuWireFramePointCloud(cloud_part_wire_frame);
//
//	// ��תƽ������Դ����
//	Eigen::Affine3f transform = Eigen::Affine3f::Identity();
//	Eigen::Vector3f ANGLE_origin;
//	Eigen::Vector3f TRANSLATE_origin;
//	ANGLE_origin << M_PI / 6, M_PI / 8, -M_PI / 3;
//	TRANSLATE_origin << 50.0, -30.0, 70.0;
//	transform.translation() << TRANSLATE_origin;
//	transform.rotate(Eigen::AngleAxisf(ANGLE_origin.z(), Eigen::Vector3f::UnitZ()));
//	transform.rotate(Eigen::AngleAxisf(ANGLE_origin.x(), Eigen::Vector3f::UnitX()));
//	transform.rotate(Eigen::AngleAxisf(ANGLE_origin.y(), Eigen::Vector3f::UnitY()));
//	pcl::transformPointCloud(*cloud_part_wire_frame, *cloud_part_wire_frame, transform);
//
//	// ������淨��
//	clock_t fpfh_start_time = clock();
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_src; // Դ
//	ne_src.setInputCloud(cloud_part_wire_frame);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_src.setSearchMethod(tree_src);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_src.setRadiusSearch(1.8);
//	ne_src.setNumberOfThreads(8);
//	ne_src.compute(*cloud_src_normals);
//
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_tar; // Ŀ��
//	ne_tar.setInputCloud(cloud_wire_frame);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_tar(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_tar.setSearchMethod(tree_tar);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_tar_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_tar.setRadiusSearch(1.8);
//	ne_src.setNumberOfThreads(8);
//	ne_tar.compute(*cloud_tar_normals);
//
//	//����FPFH
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_src;
//	fpfh_src.setInputCloud(cloud_part_wire_frame);
//	fpfh_src.setInputNormals(cloud_src_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_src_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_src.setSearchMethod(tree_src_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_src(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_src.setRadiusSearch(1.6);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_src.compute(*fpfhs_src);
//	std::cout << "compute *cloud_src fpfh" << endl;
//
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_tar;
//	fpfh_tar.setInputCloud(cloud_wire_frame);
//	fpfh_tar.setInputNormals(cloud_tar_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_tgt_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_tar.setSearchMethod(tree_tgt_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_tar(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_tar.setRadiusSearch(1.6);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_tar.compute(*fpfhs_tar);
//
//	std::cout << "compute *cloud_tar fpfh" << endl;
//	clock_t fpfh_end_time = clock();
//	cout << "fpfh time: " << (double)(fpfh_end_time - fpfh_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//SAC��׼
//	clock_t sac_start_time = clock();
//
//	pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> scia;
//	scia.setInputSource(cloud_part_wire_frame);
//	scia.setInputTarget(cloud_wire_frame);
//	scia.setSourceFeatures(fpfhs_src);
//	scia.setTargetFeatures(fpfhs_tar);
//	pcl::PointCloud<pcl::PointXYZ>::Ptr sac_result(new pcl::PointCloud<pcl::PointXYZ>);
//	scia.align(*sac_result);
//	Eigen::Matrix4f sac_trans = scia.getFinalTransformation();
//	clock_t sac_end_time = clock();
//	cout << "sac time: " << (double)(sac_end_time - sac_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//icp��׼
//	clock_t icp_start_time = clock();
//	pcl::PointCloud<pcl::PointXYZ>::Ptr icp_result(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
//	icp.setInputSource(cloud_part_wire_frame);
//	icp.setInputTarget(cloud_wire_frame);
//	// ��������Ӧ����� (e.g., correspondences with higher distances will be ignored)
//	icp.setMaxCorrespondenceDistance(40);
//	// ����������
//	icp.setMaximumIterations(100);
//	// ���α仯����֮��Ĳ�ֵ
//	icp.setTransformationEpsilon(1e-10);
//	// �������
//	icp.setEuclideanFitnessEpsilon(0.2);
//	icp.align(*icp_result, sac_trans);
//	clock_t icp_end_time = clock();
//	cout << "icp time: " << (double)(icp_end_time - icp_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	// ��ʾ����
//	// ��ɫ��Դ����
//	// ��ɫ��Ŀ�����
//	// ��ɫ����׼����
//	pcl::visualization::PCLVisualizer viewer("registration Viewer");
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h(cloud_part_wire_frame, 0, 255, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> tgt_h(cloud_wire_frame, 255, 0, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> final_h(icp_result, 0, 0, 255);
//	viewer.setBackgroundColor(0.15, 0.15, 0.15);
//
//	viewer.addPointCloud(cloud_part_wire_frame, src_h, "tar_points");
//	viewer.addPointCloud(cloud_wire_frame, tgt_h, "src_points");
//	viewer.addPointCloud(icp_result, final_h, "icp_points");
//
//	viewer.addCoordinateSystem(20.0); // �������ָʾ�ᣬ���ű���
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(100);
//		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	}
//	return 0;
//
//}
//
//// ����T����ĸ�߿����
//void createTluomuWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame)
//{
//
//	// �����µ׺�����
//	for (int i = 0; i <= 40; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i; p1.y = 0; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = i; p2.y = 50; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//	}
//
//	// ���ɶ�����
//	for (int i = 0; i <= 10; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i; p1.y = 0; p1.z = 15;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 30 + i; p2.y = 0; p2.z = 15;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = i; p3.y = 50; p3.z = 15;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 30 + i; p4.y = 50; p4.z = 15;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//	// �����ϲ������
//	for (int i = 0; i <= 20; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i + 10; p1.y = 0; p1.z = 30;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = i + 10; p2.y = 50; p2.z = 30;
//		cloud_wire_frame->points.push_back(p2);
//	}
//
//
//	// ���������ݱ�
//	for (int i = 0; i <= 50; i++)
//	{
//		// �ײ�
//		pcl::PointXYZ p1;
//		p1.x = 0; p1.y = i; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 40; p2.y = i; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//
//		// ����
//		pcl::PointXYZ p3;
//		p3.x = 0; p3.y = i; p3.z = 15;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 40; p4.y = i; p4.z = 15;
//		cloud_wire_frame->points.push_back(p4);
//
//		pcl::PointXYZ p5;
//		p5.x = 10; p5.y = i; p5.z = 15;
//		cloud_wire_frame->points.push_back(p5);
//
//		pcl::PointXYZ p6;
//		p6.x = 30; p6.y = i; p6.z = 15;
//		cloud_wire_frame->points.push_back(p6);
//
//		// �ϲ�
//		pcl::PointXYZ p7;
//		p7.x = 10; p7.y = i; p7.z = 30;
//		cloud_wire_frame->points.push_back(p7);
//
//		pcl::PointXYZ p8;
//		p8.x = 30; p8.y = i; p8.z = 30;
//		cloud_wire_frame->points.push_back(p8);
//
//	}
//
//	// ���ɵײ�������
//	for (int i = 0; i <= 15; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 0; p1.y = 0; p1.z = i;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 40; p2.y = 0; p2.z = i;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 40; p3.y = 50; p3.z = i;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 0; p4.y = 50; p4.z = i;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//	// �����ϲ�������
//	for (int i = 0; i <= 15; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 10; p1.y = 0; p1.z = i + 15;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 30; p2.y = 0; p2.z = i + 15;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 10; p3.y = 50; p3.z = i + 15;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 30; p4.y = 50; p4.z = i + 15;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//}
//
//// ���ɲ���T����ĸ�߿����
//void createPartTluomuWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame)
//{
//
//	// �����µ׺�����
//	for (int i = 0; i <= 40; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i; p1.y = 0; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = i; p2.y = 50; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//	}
//
//	// ���ɶ�����
//	for (int i = 0; i <= 10; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i; p1.y = 0; p1.z = 15;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 30 + i; p2.y = 0; p2.z = 15;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = i; p3.y = 50; p3.z = 15;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 30 + i; p4.y = 50; p4.z = 15;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//	// �����ϲ������
//	for (int i = 0; i <= 20; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i + 10; p1.y = 0; p1.z = 30;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = i + 10; p2.y = 50; p2.z = 30;
//		cloud_wire_frame->points.push_back(p2);
//	}
//
//
//	// ���������ݱ�
//	for (int i = 0; i <= 50; i++)
//	{
//		// �ײ�
//		pcl::PointXYZ p1;
//		p1.x = 0; p1.y = i; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 40; p2.y = i; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//
//		// ����
//		pcl::PointXYZ p3;
//		p3.x = 0; p3.y = i; p3.z = 15;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 40; p4.y = i; p4.z = 15;
//		cloud_wire_frame->points.push_back(p4);
//
//		pcl::PointXYZ p5;
//		p5.x = 10; p5.y = i; p5.z = 15;
//		cloud_wire_frame->points.push_back(p5);
//
//		pcl::PointXYZ p6;
//		p6.x = 30; p6.y = i; p6.z = 15;
//		cloud_wire_frame->points.push_back(p6);
//
//		// �ϲ�
//		pcl::PointXYZ p7;
//		p7.x = 10; p7.y = i; p7.z = 30;
//		cloud_wire_frame->points.push_back(p7);
//
//		pcl::PointXYZ p8;
//		p8.x = 30; p8.y = i; p8.z = 30;
//		cloud_wire_frame->points.push_back(p8);
//
//	}
//
//	// ���ɵײ�������
//	for (int i = 0; i <= 15; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 0; p1.y = 0; p1.z = i;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 40; p2.y = 0; p2.z = i;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 40; p3.y = 50; p3.z = i;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 0; p4.y = 50; p4.z = i;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//	// �����ϲ�������
//	for (int i = 0; i <= 15; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 10; p1.y = 0; p1.z = i + 15;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 30; p2.y = 0; p2.z = i + 15;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 10; p3.y = 50; p3.z = i + 15;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 30; p4.y = 50; p4.z = i + 15;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//}



///* �������߿������׼ */
//#include "TestC.h"
//#include <iostream>
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/common/io.h>
//#include <pcl/point_cloud.h>
//
//#include <pcl/PolygonMesh.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <pcl/io/vtk_lib_io.h>
//#include <pcl/io/io.h>
//#include <pcl/io/pcd_io.h>//pcd ��д����ص�ͷ�ļ���
//#include <pcl/point_types.h> //PCL��֧�ֵĵ�����ͷ�ļ�
//#include <pcl/registration/ia_ransac.h>
//#include <pcl/point_cloud.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/features/normal_3d_omp.h>
//#include <pcl/features/fpfh.h>
//#include <pcl/features/fpfh_omp.h>
//#include <pcl/search/kdtree.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/filter.h>
//#include <pcl/registration/icp.h>
//#include <pcl/visualization/pcl_visualizer.h>
//#include <time.h>
//
//int user_data;
//using namespace std;
//
//void createCuboidWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame);
//void createPartCuboidWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame);
//
//int main()
//{
//
//	/* �߿���߿������׼ */
//
//		// ������ʾ����
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame(new pcl::PointCloud<pcl::PointXYZ>());
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_part_wire_frame(new pcl::PointCloud<pcl::PointXYZ>());
//	createCuboidWireFramePointCloud(cloud_wire_frame);
//	createPartCuboidWireFramePointCloud(cloud_part_wire_frame);
//
//	// ��תƽ������Դ����
//	Eigen::Affine3f transform = Eigen::Affine3f::Identity();
//	Eigen::Vector3f ANGLE_origin;
//	Eigen::Vector3f TRANSLATE_origin;
//	ANGLE_origin << M_PI / 6, M_PI / 8, -M_PI / 3;
//	TRANSLATE_origin << 50.0, -30.0, 70.0;
//	transform.translation() << TRANSLATE_origin;
//	transform.rotate(Eigen::AngleAxisf(ANGLE_origin.z(), Eigen::Vector3f::UnitZ()));
//	transform.rotate(Eigen::AngleAxisf(ANGLE_origin.x(), Eigen::Vector3f::UnitX()));
//	transform.rotate(Eigen::AngleAxisf(ANGLE_origin.y(), Eigen::Vector3f::UnitY()));
//	pcl::transformPointCloud(*cloud_part_wire_frame, *cloud_part_wire_frame, transform);
//
//	// ������淨��
//	clock_t fpfh_start_time = clock();
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_src; // Դ
//	ne_src.setInputCloud(cloud_part_wire_frame);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_src.setSearchMethod(tree_src);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_src.setRadiusSearch(1.8);
//	ne_src.setNumberOfThreads(8);
//	ne_src.compute(*cloud_src_normals);
//
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_tar; // Ŀ��
//	ne_tar.setInputCloud(cloud_wire_frame);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_tar(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_tar.setSearchMethod(tree_tar);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_tar_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_tar.setRadiusSearch(1.8);
//	ne_src.setNumberOfThreads(8);
//	ne_tar.compute(*cloud_tar_normals);
//
//	//����FPFH
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_src;
//	fpfh_src.setInputCloud(cloud_part_wire_frame);
//	fpfh_src.setInputNormals(cloud_src_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_src_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_src.setSearchMethod(tree_src_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_src(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_src.setRadiusSearch(1.6);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_src.compute(*fpfhs_src);
//	std::cout << "compute *cloud_src fpfh" << endl;
//
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_tar;
//	fpfh_tar.setInputCloud(cloud_wire_frame);
//	fpfh_tar.setInputNormals(cloud_tar_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_tgt_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_tar.setSearchMethod(tree_tgt_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_tar(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_tar.setRadiusSearch(1.6);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_tar.compute(*fpfhs_tar);
//
//	std::cout << "compute *cloud_tar fpfh" << endl;
//	clock_t fpfh_end_time = clock();
//	cout << "fpfh time: " << (double)(fpfh_end_time - fpfh_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//SAC��׼
//	clock_t sac_start_time = clock();
//
//	pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> scia;
//	scia.setInputSource(cloud_part_wire_frame);
//	scia.setInputTarget(cloud_wire_frame);
//	scia.setSourceFeatures(fpfhs_src);
//	scia.setTargetFeatures(fpfhs_tar);
//
//	//scia.setMinSampleDistance(1);
//	//scia.setNumberOfSamples(2);
//	//scia.setCorrespondenceRandomness(20);
//
//	pcl::PointCloud<pcl::PointXYZ>::Ptr sac_result(new pcl::PointCloud<pcl::PointXYZ>);
//	scia.align(*sac_result);
//	Eigen::Matrix4f sac_trans = scia.getFinalTransformation();
//	clock_t sac_end_time = clock();
//	cout << "sac time: " << (double)(sac_end_time - sac_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//icp��׼
//	clock_t icp_start_time = clock();
//	pcl::PointCloud<pcl::PointXYZ>::Ptr icp_result(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
//	icp.setInputSource(cloud_part_wire_frame);
//	icp.setInputTarget(cloud_wire_frame);
//	// ��������Ӧ����� (e.g., correspondences with higher distances will be ignored)
//	icp.setMaxCorrespondenceDistance(40);
//	// ����������
//	icp.setMaximumIterations(100);
//	// ���α仯����֮��Ĳ�ֵ
//	icp.setTransformationEpsilon(1e-10);
//	// �������
//	icp.setEuclideanFitnessEpsilon(0.2);
//	icp.align(*icp_result, sac_trans);
//	clock_t icp_end_time = clock();
//	cout << "icp time: " << (double)(icp_end_time - icp_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	// ��ʾ����
//	// ��ɫ��Դ����
//	// ��ɫ��Ŀ�����
//	// ��ɫ����׼����
//	pcl::visualization::PCLVisualizer viewer("registration Viewer");
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h(cloud_part_wire_frame, 0, 255, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> tgt_h(cloud_wire_frame, 255, 0, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> final_h(icp_result, 0, 0, 255);
//	viewer.setBackgroundColor(0.15, 0.15, 0.15);
//
//	viewer.addPointCloud(cloud_part_wire_frame, src_h, "tar_points");
//	viewer.addPointCloud(cloud_wire_frame, tgt_h, "src_points");
//	viewer.addPointCloud(icp_result, final_h, "icp_points");
//
//
//	//�������
//	cloud_part_wire_frame->width = 1;
//	cloud_part_wire_frame->height = cloud_part_wire_frame->points.size();
//	cloud_wire_frame->width = 1;
//	cloud_wire_frame->height = cloud_wire_frame->points.size();
//	icp_result->width = 1;
//	icp_result->height = icp_result->points.size();
//	pcl::io::savePCDFileASCII("cloud_part_wire_frame.pcd", *cloud_part_wire_frame);
//	pcl::io::savePCDFileASCII("cloud_wire_frame.pcd", *cloud_wire_frame);
//	pcl::io::savePCDFileASCII("icp_result.pcd", *icp_result);
//	cout << "save point cloud file success..." << endl;
//
//	viewer.addCoordinateSystem(20.0); // �������ָʾ�ᣬ���ű���
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(100);
//		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	}
//	return 0;
//
//}
//
//
//
//// �����������߿����
//void createCuboidWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame)
//{
//
//	// ���ɳ�
//	for (int i = 0; i <= 50; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i; p1.y = 0; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = i; p2.y = 40; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = i; p3.y = 0; p3.z = 30;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = i; p4.y = 40; p4.z = 30;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//
//	// ���ɿ�
//	for (int i = 0; i <= 40; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 0; p1.y = i; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 50; p2.y = i; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 0; p3.y = i; p3.z = 30;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 50; p4.y = i; p4.z = 30;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//	// ���ɸ�
//	for (int i = 0; i <= 30; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 0; p1.y = 0; p1.z = i;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 50; p2.y = 0; p2.z = i;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 0; p3.y = 40; p3.z = i;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 50; p4.y = 40; p4.z = i;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//}
//
//// ���ɲ����������߿����
//void createPartCuboidWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame)
//{
//
//	// ���ɳ�
//	for (int i = 0; i <= 50; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i; p1.y = 0; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = i; p2.y = 40; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = i; p3.y = 0; p3.z = 30;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = i; p4.y = 40; p4.z = 30;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//
//	// ���ɿ�
//	for (int i = 0; i <= 40; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 0; p1.y = i; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 50; p2.y = i; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 0; p3.y = i; p3.z = 30;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 50; p4.y = i; p4.z = 30;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//	// ���ɸ�
//	for (int i = 0; i <= 30; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 0; p1.y = 0; p1.z = i;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 50; p2.y = 0; p2.z = i;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 0; p3.y = 40; p3.z = i;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 50; p4.y = 40; p4.z = i;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//}




///*
//  ���ɳ���
//*/
//
//#pragma once
//#include "LineAlgorithm.h"
//
//#include <iostream>  
//#include <time.h>
//
//#include <opencv2/opencv.hpp>
//#include <opencv2/ximgproc.hpp>
//#include <opencv2/xfeatures2d.hpp>
//
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/common/io.h>
//#include <pcl/point_cloud.h>
//
//#include <pcl/PolygonMesh.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <pcl/io/vtk_lib_io.h>
//#include <pcl/io/io.h>
//#include <pcl/io/pcd_io.h>//pcd ��д����ص�ͷ�ļ���
//#include <pcl/point_types.h> //PCL��֧�ֵĵ�����ͷ�ļ�
//#include <pcl/registration/ia_ransac.h>
//#include <pcl/point_cloud.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/features/normal_3d_omp.h>
//#include <pcl/features/fpfh.h>
//#include <pcl/features/fpfh_omp.h>
//#include <pcl/search/kdtree.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/filter.h>
//#include <pcl/registration/icp.h>
//#include <pcl/visualization/pcl_visualizer.h>
//
//using namespace cv;
//using namespace cv::ximgproc;
//using namespace std;
//
//// ͼƬ��С
//const int img_Width = 1280;
//const int img_Height = 720;
//
////У����ת����R��ͶӰ����P ��ͶӰ����Q
//Mat Rl, Rr, Pl, Pr, Q;
//
//// ����ģ�����
//void createCuboidWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame);
//
//// ����ת��תƽ������
//void matrix2Angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle);
//void matrix2Translate(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_translate);
//
//int main()
//{
//
//	Mat left_and_right = imread("./left_right_img/IMG_0003.jpg");
//	Rect rect_left(0, 0, img_Width, img_Height);
//	Rect rect_right(img_Width, 0, img_Width, img_Height);
//
//	Mat left = left_and_right(rect_left);
//	Mat right = left_and_right(rect_right);
//
//	//imshow("src_left", left);
//	//imshow("src_right", right);
//	//waitKey(0);
//
//	//Mat left = imread("steroMatchImage/left18.jpg");
//	//Mat right = imread("steroMatchImage/right18.jpg");
//
//	Mat cameraMatrix[2], distCoeffs[2];
//	Mat Q;
//
//	// ����ڲ�
//	cameraMatrix[0] = Mat((Mat_<double>(3, 3) << 731.5022, 0.0, 0.0,
//		0.0, 732.0644, 0.0,
//		674.6773, 357.5894, 1.0)).t();
//
//	cameraMatrix[1] = Mat((Mat_<double>(3, 3) << 727.3993, 0.0, 0.0,
//		0.0, 727.6805, 0.0,
//		634.0939, 339.6110, 1.0)).t();
//
//	// ����ϵ��
//	distCoeffs[0] = (Mat_<double>(4, 1) << 0.1070, -0.1208, 0.0, -0.0004);
//	distCoeffs[1] = (Mat_<double>(4, 1) << 0.1145, -0.1517, 0.0007, -0.0004);
//
//	// ������λ��
//	Mat R = Mat((Mat_<double>(3, 3) << 1.0, 4.7372e-04, 4.8307e-04,
//		-4.7146e-04, 1.0, -0.0047,
//		-0.0005, 0.0047, 1.0)).t();
//	Mat T = (Mat_<double>(3, 1) << -60.1450, -0.0701, -0.0956);
//
//	// ͼ��У��֮�����Ч����
//	Rect validRoi[2];
//
//	// bougust����У��
//	stereoRectify(cameraMatrix[0], distCoeffs[0],
//		cameraMatrix[1], distCoeffs[1],
//		Size(img_Width, img_Height), R, T, Rl, Rr, Pl, Pr, Q,
//		CALIB_ZERO_DISPARITY, 1, Size(img_Width, img_Height), &validRoi[0], &validRoi[1]); //validPixROI1  validPixROI2 - ��ѡ�����������Rect�����ݡ����ڲ����������ض���Ч
//
//	//��ʾ���������ͶӰ����
//	cout << "camL" << cameraMatrix[0] << endl;
//	cout << "camR" << cameraMatrix[1] << endl;
//	cout << "Pro1" << Pl << endl;
//	cout << "Pro2" << Pr << endl;
//	cout << "Q" << Q << endl;
//
//	// ����������ͼ��У������ӳ���
//	Mat rmap[2][2];
//	initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], Rl, Pl, Size(img_Width, img_Height), CV_16SC2, rmap[0][0], rmap[0][1]);
//	initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], Rr, Pr, Size(img_Width, img_Height), CV_16SC2, rmap[1][0], rmap[1][1]);
//
//	// ��Դͼ�������ӳ�䣬ӳ�����ΪУ�����ͼ��
//	Mat rleftimg, rrightimg, rcolor;
//	remap(left, rleftimg, rmap[0][0], rmap[0][1], INTER_LINEAR);
//	remap(right, rrightimg, rmap[1][0], rmap[1][1], INTER_LINEAR);
//
//	////// ����У����ͼƬ
//	////imwrite("rleft.jpg", rleftimg);
//	////imwrite("rright.jpg", rrightimg);
//
//	//// �������߽��жԱ�
//	//Mat canvas;
//	//canvas.create(img_Height, img_Width * 2, CV_8UC3);
//	//hconcat(rleftimg, rrightimg, canvas);
//	//for (int j = 0; j < canvas.rows; j += 16)
//	//	line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
//	//line(canvas, Point(1280, 0), Point(1280, 720), Scalar(255, 255, 0), 1, 8);
//
//	//// ��ʾУ�����
//	//resize(canvas, canvas, Size(), 0.5, 0.5);
//	//namedWindow("rectifyResult", WINDOW_AUTOSIZE);
//	//imshow("rectifyResult", canvas);
//	//waitKey(0);
//
//
//	// �ҶȻ�
//	Mat left_img_gray, right_img_gray;
//	cvtColor(rleftimg, left_img_gray, COLOR_BGR2GRAY);	
//	cvtColor(rrightimg, right_img_gray, COLOR_BGR2GRAY);
//	// ��
//	Mat left_img_sharpen, right_img_sharpen;
//	Mat kernel(3, 3, CV_32F, cv::Scalar(0));
//	kernel.at<float>(1, 1) = 5.0;
//	kernel.at<float>(0, 1) = -1.0;
//	kernel.at<float>(1, 0) = -1.0;
//	kernel.at<float>(1, 2) = -1.0;
//	kernel.at<float>(2, 1) = -1.0;
//	filter2D(left_img_gray, left_img_sharpen, left_img_gray.depth(), kernel);
//	filter2D(right_img_gray, right_img_sharpen, right_img_gray.depth(), kernel);
//
//	// YOLOv3 Ŀ����
//	Rect rect_left_tar(575, 138, 302, 282);
//	Rect rect_right_tar(320, 120, 306, 281);
//	Mat left_tar_img = left_img_sharpen(rect_left_tar); // ��ȡĿ������
//	Mat right_tar_img = right_img_sharpen(rect_right_tar);
//
//	//// ����Ŀ������
//	//imshow("left_lines", left_tar_img);
//	//imshow("right_lines", right_tar_img);
//	//imwrite("roi1.jpg", left_tar_img);
//	//imwrite("roi2.jpg", right_tar_img);
//	//waitKey(0);
//
//	// ��Ŀ���������ֱ�߼��
//	int    length_threshold = 30; // ���ڴ���ֵ��ֱ�߱�����
//	float  distance_threshold = 1.41421356f;  // ������߶θ�Զ�ĵ㱻��Ϊ����Ⱥ��
//	double canny_th1 = 40.0; // �ݶ�С�ڴ���ֵ����Ϊ�Ǳ�Ե
//	double canny_th2 = 100.0; // �ݶȴ��ڴ���ֵһ���Ǳ�Ե
//	int    canny_aperture_size = 3;
//	bool   do_merge = true;
//	Ptr<FastLineDetector> fld = createFastLineDetector(
//		length_threshold,
//		distance_threshold,
//		canny_th1,
//		canny_th2,
//		canny_aperture_size,
//		do_merge);
//
//	vector<Vec4f> left_lines_fld;
//	vector<Vec4f> right_lines_fld;
//	fld->detect(left_tar_img, left_lines_fld);
//	fld->detect(right_tar_img, right_lines_fld);
//
//	//// ����fld�����
//	//fld->drawSegments(left_tar_img, left_lines_fld, false);
//	//fld->drawSegments(right_tar_img, right_lines_fld, false);
//	//imshow("left_lines", left_tar_img);
//	//imshow("right_lines", right_tar_img);
//	//waitKey(0);
//
//
//	// ֱ�߸�ʽת��
//	vector<Line2D> vec_left_tar_line2d;
//	vector<Line2D> vec_right_tar_line2d;
//	for (const auto& line_4f : left_lines_fld)
//	{
//		static int serial_number = 0;
//		vec_left_tar_line2d.push_back(Line2D(serial_number, line_4f));
//		++serial_number;
//	}
//	for (const auto& line_4f : right_lines_fld)
//	{
//		static int serial_number = 0;
//		vec_right_tar_line2d.push_back(Line2D(serial_number, line_4f));
//		++serial_number;
//	}
//	vec_left_tar_line2d.push_back(Line2D(8, Vec4f(145, 188, 135, 250))); // push��ⲻ�����ĵ�
//	vec_right_tar_line2d.push_back(Line2D(8, Vec4f(160, 207, 181, 267)));
//
//	//// ����ֱ�ߺ�������
//	//Mat left_line_img = left_tar_img.clone();
//	//Mat right_line_img = right_tar_img.clone();
//	//cvtColor(left_line_img, left_line_img, COLOR_GRAY2BGR);
//	//cvtColor(right_line_img, right_line_img, COLOR_GRAY2BGR);
//	//for (const auto& line2d : vec_left_tar_line2d)
//	//{
//	//	cv::line(left_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(left_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_tar_line2d)
//	//{
//	//	cv::line(right_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(right_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//imshow("left_lines", left_line_img);
//	//imshow("right_lines", right_line_img);
//	//waitKey(0);
//
//
//	// ��ֱ�߶˵�ȡ����
//	vector<int> left_line_points_index;
//	vector<int> right_line_points_index;
//	vector<LinePoint2d> vec_left_line_points; // ������ֱ�߶˵�
//	vector<LinePoint2d> vec_right_line_points; // ������ֱ�߶˵�
//	for (const auto& line2d : vec_left_tar_line2d) // ��
//	{
//		static int index = 0;
//
//		LinePoint2d line_point1(line2d.start_point, index, line2d.serial_number, true);
//		LinePoint2d line_point2(line2d.end_point, index + 1, line2d.serial_number, false);
//
//		vec_left_line_points.push_back(line_point1);
//		vec_left_line_points.push_back(line_point2);
//		left_line_points_index.push_back(index);
//		left_line_points_index.push_back(index + 1);
//
//		index += 2;
//	}
//	for (const auto& line2d : vec_right_tar_line2d) // ��
//	{
//		static int index = 0;
//
//		LinePoint2d line_point1(line2d.start_point, index, line2d.serial_number, true);
//		LinePoint2d line_point2(line2d.end_point, index + 1, line2d.serial_number, false);
//
//		vec_right_line_points.push_back(line_point1);
//		vec_right_line_points.push_back(line_point2);
//		right_line_points_index.push_back(index);
//		right_line_points_index.push_back(index + 1);
//
//		index += 2;
//	}
//
//	// �Զ˵������ֱ�߽��з���
//	vector<vector<int>> vec_left_near_points_group; // �洢����vec_left_line_points���±�
//	vector<vector<int>> vec_right_near_points_group;
//
//	double line_point_near_threshold = 25; // �˵��������ֵ
//
//	while (left_line_points_index.size() != 0)
//	{
//		int curr_index = left_line_points_index[0];
//		vector<int> curr_near;
//		LinePoint2d curr_point = vec_left_line_points[curr_index];
//		for (int j = 0; j < left_line_points_index.size(); j++)
//		{
//			int comp_index = left_line_points_index[j];
//
//			LinePoint2d comp_point = vec_left_line_points[comp_index];
//
//			double dis = computeTowPointDistance(curr_point.point_2d, comp_point.point_2d);
//
//			if (dis <= line_point_near_threshold)
//			{
//				curr_near.push_back(comp_index);
//			}
//		}
//
//		// ɾ���Ѿ���ԵĶ˵�
//		for (int i = 0; i < curr_near.size(); i++)
//		{
//			deleteElement(left_line_points_index, curr_near[i]);
//		}
//		vec_left_near_points_group.push_back(curr_near);
//	}
//
//	while (right_line_points_index.size() != 0)
//	{
//		int curr_index = right_line_points_index[0];
//		vector<int> curr_near;
//		LinePoint2d curr_point = vec_right_line_points[curr_index];
//		for (int j = 0; j < right_line_points_index.size(); j++)
//		{
//			int comp_index = right_line_points_index[j];
//
//			LinePoint2d comp_point = vec_right_line_points[comp_index];
//
//			double dis = computeTowPointDistance(curr_point.point_2d, comp_point.point_2d);
//
//			if (dis <= line_point_near_threshold)
//			{
//				curr_near.push_back(comp_index);
//			}
//		}
//
//		// ɾ���Ѿ���ԵĶ˵�
//		for (int i = 0; i < curr_near.size(); i++)
//		{
//			deleteElement(right_line_points_index, curr_near[i]);
//		}
//		vec_right_near_points_group.push_back(curr_near);
//	}
//
//	// ֱ�ߺϲ�
//	double line_slope_similarity_threshold = 0.5; // ֱ��б�����ƶ���ֵ
//
//	// �Զ˵������ֱ�߽����ӳ����󽻵�
//	Mat left_cross_point_image = left_tar_img.clone();
//	Mat right_cross_point_image = right_tar_img.clone();
//	cvtColor(left_cross_point_image, left_cross_point_image, COLOR_GRAY2BGR);
//	cvtColor(right_cross_point_image, right_cross_point_image, COLOR_GRAY2BGR);
//	for (const auto& near_points : vec_left_near_points_group) // ��
//	{
//		if (near_points.size() == 1)
//			continue;
//		if (near_points.size() == 2)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			LinePoint2d point1 = vec_left_line_points[point_index1];
//			LinePoint2d point2 = vec_left_line_points[point_index2];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//		    Line2D& line1 = vec_left_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_left_tar_line2d[line_serial_number2];
//			Point2f cross_point = conputeIntersectionPoint(line1, line2);
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			cv::circle(left_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//		if (near_points.size() == 3)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			int point_index3 = near_points[2];
//			LinePoint2d point1 = vec_left_line_points[point_index1];
//			LinePoint2d point2 = vec_left_line_points[point_index2];
//			LinePoint2d point3 = vec_left_line_points[point_index3];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			int line_serial_number3 = point3.line_serial_number;
//			Line2D& line1 = vec_left_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_left_tar_line2d[line_serial_number2];
//			Line2D& line3 = vec_left_tar_line2d[line_serial_number3];
//			Point2f cross_point1 = conputeIntersectionPoint(line1, line2);
//			Point2f cross_point2 = conputeIntersectionPoint(line1, line3);
//			Point2f cross_point3 = conputeIntersectionPoint(line2, line3);
//			Point2f cross_point = (cross_point1 + cross_point2 + cross_point3) / 3.0;
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			if (point3.is_start_point == true)
//				line3.start_point = cross_point;
//			else
//				line3.end_point = cross_point;
//
//			cv::circle(left_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//	}
//	for (const auto& near_points : vec_right_near_points_group) // ��
//	{
//		if (near_points.size() == 1)
//			continue;
//		if (near_points.size() == 2)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			LinePoint2d point1 = vec_right_line_points[point_index1];
//			LinePoint2d point2 = vec_right_line_points[point_index2];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			Line2D& line1 = vec_right_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_right_tar_line2d[line_serial_number2];
//			Point2f cross_point = conputeIntersectionPoint(line1, line2);
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			cv::circle(right_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//		if (near_points.size() == 3)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			int point_index3 = near_points[2];
//			LinePoint2d point1 = vec_right_line_points[point_index1];
//			LinePoint2d point2 = vec_right_line_points[point_index2];
//			LinePoint2d point3 = vec_right_line_points[point_index3];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			int line_serial_number3 = point3.line_serial_number;
//			Line2D& line1 = vec_right_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_right_tar_line2d[line_serial_number2];
//			Line2D& line3 = vec_right_tar_line2d[line_serial_number3];
//			Point2f cross_point1 = conputeIntersectionPoint(line1, line2);
//			Point2f cross_point2 = conputeIntersectionPoint(line1, line3);
//			Point2f cross_point3 = conputeIntersectionPoint(line2, line3);
//			Point2f cross_point = (cross_point1 + cross_point2 + cross_point3) / 3.0;
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			if (point3.is_start_point == true)
//				line3.start_point = cross_point;
//			else
//				line3.end_point = cross_point;
//
//			cv::circle(right_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//	}
//
//	//// ��ʾֱ�߽���
//	//imshow("left_cross_point_image", left_cross_point_image);
//	//imshow("right_cross_point_image", right_cross_point_image);
//	//waitKey(0);
//
//	//// ��ʾ����������Ϊ�˵��ֱ��
//	//Mat left_line_img = left_tar_img.clone();
//	//Mat right_line_img = right_tar_img.clone();
//	//cvtColor(left_line_img, left_line_img, COLOR_GRAY2BGR);
//	//cvtColor(right_line_img, right_line_img, COLOR_GRAY2BGR);
//	//for (const auto& line2d : vec_left_tar_line2d)
//	//{
//	//	cv::line(left_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(left_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_tar_line2d)
//	//{
//	//	cv::line(right_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(right_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//imshow("left_lines", left_line_img);
//	//imshow("right_lines", right_line_img);
//	//waitKey(0);
//
//
//	// ֱ�߾ֲ�����ӳ�䵽ȫ�����꣬���¼���б��
//	vector<Line2D> vec_left_line2d;
//	vector<Line2D> vec_right_line2d;
//	for (const auto& line2d : vec_left_tar_line2d) // ��
//	{
//		static int index = 0;
//		float x1 = line2d.start_point.x + rect_left_tar.x;
//		float y1 = line2d.start_point.y + rect_left_tar.y;
//		float x2 = line2d.end_point.x + rect_left_tar.x;
//		float y2 = line2d.end_point.y + rect_left_tar.y;
//		Vec4f line4f(x1, y1, x2, y2);
//		vec_left_line2d.push_back(Line2D(index, line4f));
//		++index;
//	}
//	for (const auto& line2d : vec_right_tar_line2d) // ��
//	{
//		static int index = 0;
//		float x1 = line2d.start_point.x + rect_right_tar.x;
//		float y1 = line2d.start_point.y + rect_right_tar.y;
//		float x2 = line2d.end_point.x + rect_right_tar.x;
//		float y2 = line2d.end_point.y + rect_right_tar.y;
//		Vec4f line4f(x1, y1, x2, y2);
//		vec_right_line2d.push_back(Line2D(index, line4f));
//		++index;
//	}
//
//	//// ���ƴ���õ�ȫ��ֱ��
//	//Mat left_global_line_img = rleftimg.clone();
//	//Mat right_global_line_img = rrightimg.clone();
//
//	//for (const auto& line2d : vec_left_line2d)
//	//{
//	//	cv::line(left_global_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 2);
//	//	cv::putText(left_global_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_line2d)
//	//{
//	//	cv::line(right_global_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 2);
//	//	cv::putText(right_global_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//
//	//Mat left_right_global_line_img; // ƴ������ͼ��
//	//left_right_global_line_img.create(img_Height, img_Width * 2, CV_8UC3);
//	//hconcat(left_global_line_img, right_global_line_img, left_right_global_line_img);
//	//imshow("left_global_line_img", left_global_line_img);
//	//imshow("right_global_line_img", right_global_line_img);
//	//imshow("left_right_global_line_img", left_right_global_line_img);
//	//waitKey(0);
//
//	// ֱ��ƥ��
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_reconstruct(new pcl::PointCloud<pcl::PointXYZ>()); // ֱ�ߵ���
//	std::vector<std::vector<int>> vec_match_index_group; // ƥ������
//	lineMatch(vec_left_line2d, vec_right_line2d, vec_match_index_group);
//
//	// ֱ���ؽ�
//	vector<Line3D> vec_line3d;
//	for (const auto vec_match_index : vec_match_index_group)
//	{
//		if (vec_match_index.size() != 2)
//			continue;
//
//		int left_index = vec_match_index[0];
//		int right_index = vec_match_index[1];
//		const Line2D& left_line = vec_left_line2d[left_index];
//		const Line2D& right_line = vec_right_line2d[right_index];
//
//		// p1��㣬p2�յ�
//		const Point2f& left_p1 = left_line.start_point;
//		const Point2f& left_p2 = left_line.end_point;
//		const Point2f& right_p1 = right_line.start_point;
//		const Point2f& right_p2 = right_line.end_point;
//
//		// ����ռ�ֱ�������˵�
//		const float f = 622.812; // ���࣬���㣬���߾���
//		const float cx = 654.193;
//		const float cy = 344.569;
//		const float Tx = 60.14;
//
//		float p1_disp = left_p1.x - right_p1.x;
//
//		float p1_z_w = f * Tx / p1_disp;
//		float p1_x_w = (left_p1.x - cx)*p1_z_w / f;
//		float p1_y_w = (left_p1.y - cy)*p1_z_w / f;
//
//		float p2_disp = left_p2.x - right_p2.x;
//		float p2_z_w = f * Tx / p2_disp;
//		float p2_x_w = (left_p2.x - cx)*p2_z_w / f;
//		float p2_y_w = (left_p2.y - cy)*p2_z_w / f;
//
//		// push�ؽ��Ŀռ�ֱ��
//		Line3D line_3d(Point3f(p1_x_w, p1_y_w, p1_z_w), Point3f(p2_x_w, p2_y_w, p2_z_w));
//		vec_line3d.push_back(line_3d);
//
//		// push �ռ�ֱ�߶˵�
//		pcl::PointXYZ p1_w(p1_x_w, p1_y_w, p1_z_w);
//		pcl::PointXYZ p2_w(p2_x_w, p2_y_w, p2_z_w);
//		cloud_reconstruct->points.push_back(p1_w);
//		cloud_reconstruct->points.push_back(p2_w);
//	}
//
//	// ȫ��ֱ�߲����ɵ���
//	float sample_dis = 1.0; // ֱ�ߴ�Ͼ���
//
//	for (const auto& line_3d : vec_line3d)
//	{
//		// ֱ�ߵȾ�����
//		float x0 = line_3d.x0_vec.x;
//		float y0 = line_3d.x0_vec.y;
//		float z0 = line_3d.x0_vec.z;
//
//		float v1 = line_3d.normal_vec.x;
//		float v2 = line_3d.normal_vec.y;
//		float v3 = line_3d.normal_vec.z;
//
//		float v_square = line_3d.normal_vec.dot(line_3d.normal_vec);
//		float deta_t = sample_dis / std::sqrt(v_square);
//
//		for (float t = line_3d.t_start; t < line_3d.t_end; t += deta_t)
//		{
//			pcl::PointXYZ point_3d;
//			point_3d.x = x0 + v1 * t;
//			point_3d.y = y0 + v2 * t;
//			point_3d.z = z0 + v3 * t;
//			cloud_reconstruct->points.push_back(point_3d);
//		}
//	}
//
//	//// ȫ��ֱ�߲����ɵ���
//	//float sample_dis = 1.0; // ֱ�ߴ�Ͼ���
//
//	//for (const auto& line_3d : vec_line3d)
//	//{
//	//	// ֱ�ߵȾ�����
//	//	Point3f deta_p = (line_3d.start_point - line_3d.end_point)/ line_3d.length;
//
//	//	for (float t = 0; t < line_3d.length; ++t)
//	//	{
//	//		static pcl::PointXYZ point_3d(line_3d.start_point.x, line_3d.start_point.y, line_3d.start_point.z);
//	//		point_3d.x = point_3d.x + deta_p.x;
//	//		point_3d.y = point_3d.y + deta_p.y;
//	//		point_3d.z = point_3d.z + deta_p.z;
//	//		cloud_reconstruct->points.push_back(point_3d);
//	//	}
//	//}
//
//	// ������׼
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src(new pcl::PointCloud<pcl::PointXYZ>());
//	createCuboidWireFramePointCloud(cloud_src);
//	// ������Ʊ��淨��
//	clock_t fpfh_start_time = clock();
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_src; // Դ
//	ne_src.setInputCloud(cloud_src);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_src.setSearchMethod(tree_src);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_src.setRadiusSearch(2.0);
//	ne_src.setNumberOfThreads(8);
//	ne_src.compute(*cloud_src_normals);
//
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_tar; // �ؽ�
//	ne_tar.setInputCloud(cloud_reconstruct);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_reconstruct(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_tar.setSearchMethod(tree_reconstruct);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_reconstruct_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_tar.setRadiusSearch(2.0);
//	ne_src.setNumberOfThreads(8);
//	ne_tar.compute(*cloud_reconstruct_normals);
//
//	//����FPFH
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_src; // Դ
//	fpfh_src.setInputCloud(cloud_src);
//	fpfh_src.setInputNormals(cloud_src_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_src_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_src.setSearchMethod(tree_src_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_src(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_src.setRadiusSearch(1.8);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_src.compute(*fpfhs_src);
//	std::cout << "compute *cloud_src fpfh" << endl;
//
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_tar; // �ؽ�
//	fpfh_tar.setInputCloud(cloud_reconstruct);
//	fpfh_tar.setInputNormals(cloud_reconstruct_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_reconstruct_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_tar.setSearchMethod(tree_reconstruct_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_reconstruct(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_tar.setRadiusSearch(1.8);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_tar.compute(*fpfhs_reconstruct);
//
//	std::cout << "compute *cloud_tar fpfh" << endl;
//	clock_t fpfh_end_time = clock();
//	cout << "fpfh time: " << (double)(fpfh_end_time - fpfh_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//SAC��׼
//	clock_t sac_start_time = clock();
//
//	pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> scia;
//	scia.setInputSource(cloud_src);
//	scia.setInputTarget(cloud_reconstruct);
//	scia.setSourceFeatures(fpfhs_src);
//	scia.setTargetFeatures(fpfhs_reconstruct);
//
//	//scia.setMinSampleDistance(1);
//	//scia.setNumberOfSamples(2);
//	//scia.setCorrespondenceRandomness(20);
//
//	pcl::PointCloud<pcl::PointXYZ>::Ptr sac_result(new pcl::PointCloud<pcl::PointXYZ>);
//	scia.align(*sac_result);
//	Eigen::Matrix4f sac_trans = scia.getFinalTransformation();
//	clock_t sac_end_time = clock();
//	cout << "sac time: " << (double)(sac_end_time - sac_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//icp��׼
//	clock_t icp_start_time = clock();
//	pcl::PointCloud<pcl::PointXYZ>::Ptr icp_result(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
//	icp.setInputSource(cloud_src);
//	icp.setInputTarget(cloud_reconstruct);
//	// ��������Ӧ����� (e.g., correspondences with higher distances will be ignored)
//	icp.setMaxCorrespondenceDistance(40);
//	// ����������
//	icp.setMaximumIterations(100);
//	// ���α仯����֮��Ĳ�ֵ
//	icp.setTransformationEpsilon(1e-10);
//	// �������
//	icp.setEuclideanFitnessEpsilon(0.2);
//	icp.align(*icp_result, sac_trans);
//	clock_t icp_end_time = clock();
//	cout << "icp time: " << (double)(icp_end_time - icp_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//	// ��ȡ���յı任����
//	Eigen::Matrix4f final_trans_pcl = icp.getFinalTransformation();
//	cv::Mat_<float> final_trans_cv(4, 4);
//	for (size_t i = 0; i < 4; i++)
//	{
//		for (size_t j = 0; j < 4; j++)
//		{
//			final_trans_cv(i, j) = final_trans_pcl(i, j);
//		}
//	}
//
//	// ��������ϵ������
//	vector<cv::Point3f> vec_axis_3d;
//	vector<cv::Point2f> vec_axis_2d;
//	vec_axis_3d.push_back(Point3f(10.0, 0.0, 0.0));
//	vec_axis_3d.push_back(Point3f(0.0, 10.0, 0.0));
//	vec_axis_3d.push_back(Point3f(0.0, 0.0, 10.0));
//	vec_axis_3d.push_back(Point3f(0.0, 0.0, 0.0));
//
//	// ��������ϵת��������
//	Eigen::Vector3f final_angle;
//	Eigen::Vector3f final_translate;
//	matrix2Angle(final_trans_pcl, final_angle);
//	matrix2Translate(final_trans_pcl, final_translate);
//	cv::Mat rVec(3, 1, cv::DataType<float>::type); // Rotation vector
//	rVec.at<float>(0) = final_angle.x();
//	rVec.at<float>(1) = final_angle.y();
//	rVec.at<float>(2) = final_angle.z();
//
//	cv::Mat tVec(3, 1, cv::DataType<float>::type); // Translation vector
//	tVec.at<float>(0) = final_translate.x();
//	tVec.at<float>(1) = final_translate.y();
//	tVec.at<float>(2) = final_translate.z();
//	
//	cv::projectPoints(vec_axis_3d, rVec, tVec, Pl.colRange(0, 3), Mat(), vec_axis_2d);
//
//	// ������ͶӰ���������
//	Mat reproject_axis_img = rleftimg.clone();
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[0], Scalar(0, 0, 255), 2);
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[1], Scalar(0, 255, 0), 2);
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[2], Scalar(255, 0, 0), 2);
//	imshow("reproject_axis_img", reproject_axis_img);
//	waitKey(0);
//
//	// ��ʾ����
//	// ��ɫ��Դ����
//	// ��ɫ��Ŀ�����
//	// ��ɫ����׼����
//	pcl::visualization::PCLVisualizer viewer("registration Viewer");
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h(cloud_src, 0, 255, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> tgt_h(cloud_reconstruct, 255, 0, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> final_h(icp_result, 0, 0, 255);
//	viewer.setBackgroundColor(1.0, 1.0, 1.0);
//
//	//viewer.addPointCloud(cloud_src, src_h, "tar_points");
//	viewer.addPointCloud(cloud_reconstruct, tgt_h, "src_points");
//	viewer.addPointCloud(icp_result, final_h, "icp_points");
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "src_points"); // ���õ��ƴ�С
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "tar_points");
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "icp_points");
//	viewer.addCoordinateSystem(20.0); // �������ָʾ�ᣬ���ű���-------��ɫX����ɫY����ɫZ
//
//	//// �������
//	//cloud_reconstruct->width = 1;
//	//cloud_reconstruct->height = cloud_reconstruct->points.size();
//	//pcl::io::savePCDFileASCII("./cloud_reconstruct/cloud_reconstruct.pcd", *cloud_reconstruct);
//	//cout << "cloud_reconstruct is saved" << endl;
//
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(100);
//		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	}
//
//	return 0;
//}
//
//
//// �����������߿����
//void createCuboidWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame)
//{
//
//	// ���ɳ�
//	for (int i = 0; i <= 50; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i; p1.y = 0; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = i; p2.y = 40; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = i; p3.y = 0; p3.z = 30;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = i; p4.y = 40; p4.z = 30;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//
//	// ���ɿ�
//	for (int i = 0; i <= 40; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 0; p1.y = i; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 50; p2.y = i; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 0; p3.y = i; p3.z = 30;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 50; p4.y = i; p4.z = 30;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//	// ���ɸ�
//	for (int i = 0; i <= 30; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 0; p1.y = 0; p1.z = i;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 50; p2.y = 0; p2.z = i;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 0; p3.y = 40; p3.z = i;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 50; p4.y = 40; p4.z = i;
//		cloud_wire_frame->points.push_back(p4);
//	}
//}
//
//
////����תƽ�ƾ��������ת�Ƕ�
//void matrix2Angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle)
//{
//	double ax, ay, az;
//	if (result_trans(2, 0) == 1 || result_trans(2, 0) == -1)
//	{
//		az = 0;
//		double dlta;
//		dlta = atan2(result_trans(0, 1), result_trans(0, 2));
//		if (result_trans(2, 0) == -1)
//		{
//			ay = M_PI / 2;
//			ax = az + dlta;
//		}
//		else
//		{
//			ay = -M_PI / 2;
//			ax = -az + dlta;
//		}
//	}
//	else
//	{
//		ay = -asin(result_trans(2, 0));
//		ax = atan2(result_trans(2, 1) / cos(ay), result_trans(2, 2) / cos(ay));
//		az = atan2(result_trans(1, 0) / cos(ay), result_trans(0, 0) / cos(ay));
//	}
//	result_angle << ax, ay, az;
//}
//
////����תƽ�ƾ������ƽ�ƾ���
//void matrix2Translate(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_translate)
//{
//	result_translate << result_trans(0, 3), result_trans(1, 3), result_trans(2, 3);
//}


///*
//  L��Ƭ
//*/
//
//#pragma once
//#include "LineAlgorithm.h"
//
//#include <iostream>  
//#include <time.h>
//
//#include <opencv2/opencv.hpp>
//#include <opencv2/ximgproc.hpp>
//#include <opencv2/xfeatures2d.hpp>
//
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/common/io.h>
//#include <pcl/point_cloud.h>
//
//#include <pcl/PolygonMesh.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <pcl/io/vtk_lib_io.h>
//#include <pcl/io/io.h>
//#include <pcl/io/pcd_io.h>//pcd ��д����ص�ͷ�ļ���
//#include <pcl/point_types.h> //PCL��֧�ֵĵ�����ͷ�ļ�
//#include <pcl/registration/ia_ransac.h>
//#include <pcl/point_cloud.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/features/normal_3d_omp.h>
//#include <pcl/features/fpfh.h>
//#include <pcl/features/fpfh_omp.h>
//#include <pcl/search/kdtree.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/filter.h>
//#include <pcl/registration/icp.h>
//#include <pcl/visualization/pcl_visualizer.h>
//
//using namespace cv;
//using namespace cv::ximgproc;
//using namespace std;
//
//// ͼƬ��С
//const int img_Width = 1280;
//const int img_Height = 720;
//
////У����ת����R��ͶӰ����P ��ͶӰ����Q
//Mat Rl, Rr, Pl, Pr, Q;
//
//// ����ģ�����
//void createLWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame);
//
//// ����ת��תƽ������
//void matrix2Angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle);
//void matrix2Translate(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_translate);
//
//int main()
//{
//
//	Mat left_and_right = imread("./left_right_img/IMG_0007.jpg");
//	Rect rect_left(0, 0, img_Width, img_Height);
//	Rect rect_right(img_Width, 0, img_Width, img_Height);
//
//	Mat left = left_and_right(rect_left);
//	Mat right = left_and_right(rect_right);
//
//	//imshow("src_left", left);
//	//imshow("src_right", right);
//	//waitKey(0);
//
//	Mat cameraMatrix[2], distCoeffs[2];
//	Mat Q;
//
//	// ����ڲ�
//	cameraMatrix[0] = Mat((Mat_<double>(3, 3) << 731.5022, 0.0, 0.0,
//		0.0, 732.0644, 0.0,
//		674.6773, 357.5894, 1.0)).t();
//
//	cameraMatrix[1] = Mat((Mat_<double>(3, 3) << 727.3993, 0.0, 0.0,
//		0.0, 727.6805, 0.0,
//		634.0939, 339.6110, 1.0)).t();
//
//	// ����ϵ��
//	distCoeffs[0] = (Mat_<double>(4, 1) << 0.1070, -0.1208, 0.0, -0.0004);
//	distCoeffs[1] = (Mat_<double>(4, 1) << 0.1145, -0.1517, 0.0007, -0.0004);
//
//	// ������λ��
//	Mat R = Mat((Mat_<double>(3, 3) << 1.0, 4.7372e-04, 4.8307e-04,
//		-4.7146e-04, 1.0, -0.0047,
//		-0.0005, 0.0047, 1.0)).t();
//	Mat T = (Mat_<double>(3, 1) << -60.1450, -0.0701, -0.0956);
//
//	// ͼ��У��֮�����Ч����
//	Rect validRoi[2];
//
//	// bougust����У��
//	stereoRectify(cameraMatrix[0], distCoeffs[0],
//		cameraMatrix[1], distCoeffs[1],
//		Size(img_Width, img_Height), R, T, Rl, Rr, Pl, Pr, Q,
//		CALIB_ZERO_DISPARITY, 1, Size(img_Width, img_Height), &validRoi[0], &validRoi[1]); //validPixROI1  validPixROI2 - ��ѡ�����������Rect�����ݡ����ڲ����������ض���Ч
//
//	//��ʾ���������ͶӰ����
//	cout << "camL" << cameraMatrix[0] << endl;
//	cout << "camR" << cameraMatrix[1] << endl;
//	cout << "Pro1" << Pl << endl;
//	cout << "Pro2" << Pr << endl;
//	cout << "Q" << Q << endl;
//
//	// ����������ͼ��У������ӳ���
//	Mat rmap[2][2];
//	initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], Rl, Pl, Size(img_Width, img_Height), CV_16SC2, rmap[0][0], rmap[0][1]);
//	initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], Rr, Pr, Size(img_Width, img_Height), CV_16SC2, rmap[1][0], rmap[1][1]);
//
//	// ��Դͼ�������ӳ�䣬ӳ�����ΪУ�����ͼ��
//	Mat rleftimg, rrightimg, rcolor;
//	remap(left, rleftimg, rmap[0][0], rmap[0][1], INTER_LINEAR);
//	remap(right, rrightimg, rmap[1][0], rmap[1][1], INTER_LINEAR);
//
//	////// ����У����ͼƬ
//	////imwrite("rleft.jpg", rleftimg);
//	////imwrite("rright.jpg", rrightimg);
//
//	//// �������߽��жԱ�
//	//Mat canvas;
//	//canvas.create(img_Height, img_Width * 2, CV_8UC3);
//	//hconcat(rleftimg, rrightimg, canvas);
//	//for (int j = 0; j < canvas.rows; j += 16)
//	//	line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
//	//line(canvas, Point(1280, 0), Point(1280, 720), Scalar(255, 255, 0), 1, 8);
//
//	//// ��ʾУ�����
//	//resize(canvas, canvas, Size(), 0.5, 0.5);
//	//namedWindow("rectifyResult", WINDOW_AUTOSIZE);
//	//imshow("rectifyResult", canvas);
//	//waitKey(0);
//
//
//	// �ҶȻ�
//	Mat left_img_gray, right_img_gray;
//	cvtColor(rleftimg, left_img_gray, COLOR_BGR2GRAY);
//	cvtColor(rrightimg, right_img_gray, COLOR_BGR2GRAY);
//	// ��
//	Mat left_img_sharpen, right_img_sharpen;
//	Mat kernel(3, 3, CV_32F, cv::Scalar(0));
//	kernel.at<float>(1, 1) = 5.0;
//	kernel.at<float>(0, 1) = -1.0;
//	kernel.at<float>(1, 0) = -1.0;
//	kernel.at<float>(1, 2) = -1.0;
//	kernel.at<float>(2, 1) = -1.0;
//	filter2D(left_img_gray, left_img_sharpen, left_img_gray.depth(), kernel);
//	filter2D(right_img_gray, right_img_sharpen, right_img_gray.depth(), kernel);
//
//	// YOLOv3 Ŀ����
//	Rect rect_left_tar(485, 300, 337, 201);
//	Rect rect_right_tar(250, 321, 337, 214);
//	Mat left_tar_img = left_img_sharpen(rect_left_tar); // ��ȡĿ������
//	Mat right_tar_img = right_img_sharpen(rect_right_tar);
//
//	// //����Ŀ������
//	//imshow("left_lines", left_tar_img);
//	//imshow("right_lines", right_tar_img);
//	////imwrite("roi1.jpg", left_tar_img);
//	////imwrite("roi2.jpg", right_tar_img);
//	//waitKey(0);
//
//	// ��Ŀ���������ֱ�߼��
//	int    length_threshold = 20; // ���ڴ���ֵ��ֱ�߱�����
//	float  distance_threshold = 1.41421356f;  // ������߶θ�Զ�ĵ㱻��Ϊ����Ⱥ��
//	double canny_th1 = 40.0; // �ݶ�С�ڴ���ֵ����Ϊ�Ǳ�Ե
//	double canny_th2 = 100.0; // �ݶȴ��ڴ���ֵһ���Ǳ�Ե
//	int    canny_aperture_size = 3;
//	bool   do_merge = true;
//	Ptr<FastLineDetector> fld = createFastLineDetector(
//		length_threshold,
//		distance_threshold,
//		canny_th1,
//		canny_th2,
//		canny_aperture_size,
//		do_merge);
//
//	vector<Vec4f> left_lines_fld;
//	vector<Vec4f> right_lines_fld;
//	fld->detect(left_tar_img, left_lines_fld);
//	fld->detect(right_tar_img, right_lines_fld);
//
//	//// ����fld�����
//	//fld->drawSegments(left_tar_img, left_lines_fld, false);
//	//fld->drawSegments(right_tar_img, right_lines_fld, false);
//	//imshow("left_lines", left_tar_img);
//	//imshow("right_lines", right_tar_img);
//	//waitKey(0);
//
//
//	// ֱ�߸�ʽת��
//	vector<Line2D> vec_left_tar_line2d;
//	vector<Line2D> vec_right_tar_line2d;
//	for (const auto& line_4f : left_lines_fld)
//	{
//		static int serial_number = 0;
//		vec_left_tar_line2d.push_back(Line2D(serial_number, line_4f));
//		++serial_number;
//	}
//	for (const auto& line_4f : right_lines_fld)
//	{
//		static int serial_number = 0;
//		vec_right_tar_line2d.push_back(Line2D(serial_number, line_4f));
//		++serial_number;
//	}
//
//	//// ����ֱ�ߺ�������
//	//Mat left_line_img = left_tar_img.clone();
//	//Mat right_line_img = right_tar_img.clone();
//	//cvtColor(left_line_img, left_line_img, COLOR_GRAY2BGR);
//	//cvtColor(right_line_img, right_line_img, COLOR_GRAY2BGR);
//	//for (const auto& line2d : vec_left_tar_line2d)
//	//{
//	//	cv::line(left_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(left_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_tar_line2d)
//	//{
//	//	cv::line(right_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(right_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//imshow("left_lines", left_line_img);
//	//imshow("right_lines", right_line_img);
//	//waitKey(0);
//
//
//	// ��ֱ�߶˵�ȡ����
//	vector<int> left_line_points_index;
//	vector<int> right_line_points_index;
//	vector<LinePoint2d> vec_left_line_points; // ������ֱ�߶˵�
//	vector<LinePoint2d> vec_right_line_points; // ������ֱ�߶˵�
//	for (const auto& line2d : vec_left_tar_line2d) // ��
//	{
//		static int index = 0;
//
//		LinePoint2d line_point1(line2d.start_point, index, line2d.serial_number, true);
//		LinePoint2d line_point2(line2d.end_point, index + 1, line2d.serial_number, false);
//
//		vec_left_line_points.push_back(line_point1);
//		vec_left_line_points.push_back(line_point2);
//		left_line_points_index.push_back(index);
//		left_line_points_index.push_back(index + 1);
//
//		index += 2;
//	}
//	for (const auto& line2d : vec_right_tar_line2d) // ��
//	{
//		static int index = 0;
//
//		LinePoint2d line_point1(line2d.start_point, index, line2d.serial_number, true);
//		LinePoint2d line_point2(line2d.end_point, index + 1, line2d.serial_number, false);
//
//		vec_right_line_points.push_back(line_point1);
//		vec_right_line_points.push_back(line_point2);
//		right_line_points_index.push_back(index);
//		right_line_points_index.push_back(index + 1);
//
//		index += 2;
//	}
//
//	// �Զ˵������ֱ�߽��з���
//	vector<vector<int>> vec_left_near_points_group; // �洢����vec_left_line_points���±�
//	vector<vector<int>> vec_right_near_points_group;
//
//	double line_point_near_threshold = 25; // �˵��������ֵ
//
//	while (left_line_points_index.size() != 0)
//	{
//		int curr_index = left_line_points_index[0];
//		vector<int> curr_near;
//		LinePoint2d curr_point = vec_left_line_points[curr_index];
//		for (int j = 0; j < left_line_points_index.size(); j++)
//		{
//			int comp_index = left_line_points_index[j];
//
//			LinePoint2d comp_point = vec_left_line_points[comp_index];
//
//			double dis = computeTowPointDistance(curr_point.point_2d, comp_point.point_2d);
//
//			if (dis <= line_point_near_threshold)
//			{
//				curr_near.push_back(comp_index);
//			}
//		}
//
//		// ɾ���Ѿ���ԵĶ˵�
//		for (int i = 0; i < curr_near.size(); i++)
//		{
//			deleteElement(left_line_points_index, curr_near[i]);
//		}
//		vec_left_near_points_group.push_back(curr_near);
//	}
//
//	while (right_line_points_index.size() != 0)
//	{
//		int curr_index = right_line_points_index[0];
//		vector<int> curr_near;
//		LinePoint2d curr_point = vec_right_line_points[curr_index];
//		for (int j = 0; j < right_line_points_index.size(); j++)
//		{
//			int comp_index = right_line_points_index[j];
//
//			LinePoint2d comp_point = vec_right_line_points[comp_index];
//
//			double dis = computeTowPointDistance(curr_point.point_2d, comp_point.point_2d);
//
//			if (dis <= line_point_near_threshold)
//			{
//				curr_near.push_back(comp_index);
//			}
//		}
//
//		// ɾ���Ѿ���ԵĶ˵�
//		for (int i = 0; i < curr_near.size(); i++)
//		{
//			deleteElement(right_line_points_index, curr_near[i]);
//		}
//		vec_right_near_points_group.push_back(curr_near);
//	}
//
//	// ֱ�ߺϲ�
//	double line_slope_similarity_threshold = 0.5; // ֱ��б�����ƶ���ֵ
//
//	// �Զ˵������ֱ�߽����ӳ����󽻵�
//	Mat left_cross_point_image = left_tar_img.clone();
//	Mat right_cross_point_image = right_tar_img.clone();
//	cvtColor(left_cross_point_image, left_cross_point_image, COLOR_GRAY2BGR);
//	cvtColor(right_cross_point_image, right_cross_point_image, COLOR_GRAY2BGR);
//	for (const auto& near_points : vec_left_near_points_group) // ��
//	{
//		if (near_points.size() == 1)
//			continue;
//		if (near_points.size() == 2)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			LinePoint2d point1 = vec_left_line_points[point_index1];
//			LinePoint2d point2 = vec_left_line_points[point_index2];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			Line2D& line1 = vec_left_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_left_tar_line2d[line_serial_number2];
//			Point2f cross_point = conputeIntersectionPoint(line1, line2);
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			cv::circle(left_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//		if (near_points.size() == 3)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			int point_index3 = near_points[2];
//			LinePoint2d point1 = vec_left_line_points[point_index1];
//			LinePoint2d point2 = vec_left_line_points[point_index2];
//			LinePoint2d point3 = vec_left_line_points[point_index3];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			int line_serial_number3 = point3.line_serial_number;
//			Line2D& line1 = vec_left_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_left_tar_line2d[line_serial_number2];
//			Line2D& line3 = vec_left_tar_line2d[line_serial_number3];
//			Point2f cross_point1 = conputeIntersectionPoint(line1, line2);
//			Point2f cross_point2 = conputeIntersectionPoint(line1, line3);
//			Point2f cross_point3 = conputeIntersectionPoint(line2, line3);
//			Point2f cross_point = (cross_point1 + cross_point2 + cross_point3) / 3.0;
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			if (point3.is_start_point == true)
//				line3.start_point = cross_point;
//			else
//				line3.end_point = cross_point;
//
//			cv::circle(left_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//	}
//	for (const auto& near_points : vec_right_near_points_group) // ��
//	{
//		if (near_points.size() == 1)
//			continue;
//		if (near_points.size() == 2)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			LinePoint2d point1 = vec_right_line_points[point_index1];
//			LinePoint2d point2 = vec_right_line_points[point_index2];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			Line2D& line1 = vec_right_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_right_tar_line2d[line_serial_number2];
//			Point2f cross_point = conputeIntersectionPoint(line1, line2);
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			cv::circle(right_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//		if (near_points.size() == 3)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			int point_index3 = near_points[2];
//			LinePoint2d point1 = vec_right_line_points[point_index1];
//			LinePoint2d point2 = vec_right_line_points[point_index2];
//			LinePoint2d point3 = vec_right_line_points[point_index3];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			int line_serial_number3 = point3.line_serial_number;
//			Line2D& line1 = vec_right_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_right_tar_line2d[line_serial_number2];
//			Line2D& line3 = vec_right_tar_line2d[line_serial_number3];
//			Point2f cross_point1 = conputeIntersectionPoint(line1, line2);
//			Point2f cross_point2 = conputeIntersectionPoint(line1, line3);
//			Point2f cross_point3 = conputeIntersectionPoint(line2, line3);
//			Point2f cross_point = (cross_point1 + cross_point2 + cross_point3) / 3.0;
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			if (point3.is_start_point == true)
//				line3.start_point = cross_point;
//			else
//				line3.end_point = cross_point;
//
//			cv::circle(right_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//	}
//
//	//// ��ʾֱ�߽���
//	//imshow("left_cross_point_image", left_cross_point_image);
//	//imshow("right_cross_point_image", right_cross_point_image);
//	//waitKey(0);
//
//	//// ��ʾ����������Ϊ�˵��ֱ��
//	//Mat left_line_img = left_tar_img.clone();
//	//Mat right_line_img = right_tar_img.clone();
//	//cvtColor(left_line_img, left_line_img, COLOR_GRAY2BGR);
//	//cvtColor(right_line_img, right_line_img, COLOR_GRAY2BGR);
//	//for (const auto& line2d : vec_left_tar_line2d)
//	//{
//	//	cv::line(left_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(left_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_tar_line2d)
//	//{
//	//	cv::line(right_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(right_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//imshow("left_lines", left_line_img);
//	//imshow("right_lines", right_line_img);
//	//waitKey(0);
//
//
//	// ֱ�߾ֲ�����ӳ�䵽ȫ�����꣬���¼���б��
//	vector<Line2D> vec_left_line2d;
//	vector<Line2D> vec_right_line2d;
//	for (const auto& line2d : vec_left_tar_line2d) // ��
//	{
//		static int index = 0;
//		float x1 = line2d.start_point.x + rect_left_tar.x;
//		float y1 = line2d.start_point.y + rect_left_tar.y;
//		float x2 = line2d.end_point.x + rect_left_tar.x;
//		float y2 = line2d.end_point.y + rect_left_tar.y;
//		Vec4f line4f(x1, y1, x2, y2);
//		vec_left_line2d.push_back(Line2D(index, line4f));
//		++index;
//	}
//	for (const auto& line2d : vec_right_tar_line2d) // ��
//	{
//		static int index = 0;
//		float x1 = line2d.start_point.x + rect_right_tar.x;
//		float y1 = line2d.start_point.y + rect_right_tar.y;
//		float x2 = line2d.end_point.x + rect_right_tar.x;
//		float y2 = line2d.end_point.y + rect_right_tar.y;
//		Vec4f line4f(x1, y1, x2, y2);
//		vec_right_line2d.push_back(Line2D(index, line4f));
//		++index;
//	}
//
//	//// ���ƴ���õ�ȫ��ֱ��
//	//Mat left_global_line_img = rleftimg.clone();
//	//Mat right_global_line_img = rrightimg.clone();
//
//	//for (const auto& line2d : vec_left_line2d)
//	//{
//	//	cv::line(left_global_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 2);
//	//	cv::putText(left_global_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_line2d)
//	//{
//	//	cv::line(right_global_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 2);
//	//	cv::putText(right_global_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//
//	//Mat left_right_global_line_img; // ƴ������ͼ��
//	//left_right_global_line_img.create(img_Height, img_Width * 2, CV_8UC3);
//	//hconcat(left_global_line_img, right_global_line_img, left_right_global_line_img);
//	//imshow("left_global_line_img", left_global_line_img);
//	//imshow("right_global_line_img", right_global_line_img);
//	//imshow("left_right_global_line_img", left_right_global_line_img);
//	//waitKey(0);
//
//	// ֱ��ƥ��
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_reconstruct(new pcl::PointCloud<pcl::PointXYZ>()); // ֱ�ߵ���
//	std::vector<std::vector<int>> vec_match_index_group; // ƥ������
//	lineMatch(vec_left_line2d, vec_right_line2d, vec_match_index_group);
//
//	// ֱ���ؽ�
//	vector<Line3D> vec_line3d;
//	for (const auto vec_match_index : vec_match_index_group)
//	{
//		if (vec_match_index.size() != 2)
//			continue;
//
//		int left_index = vec_match_index[0];
//		int right_index = vec_match_index[1];
//		const Line2D& left_line = vec_left_line2d[left_index];
//		const Line2D& right_line = vec_right_line2d[right_index];
//
//		// p1��㣬p2�յ�
//		const Point2f& left_p1 = left_line.start_point;
//		const Point2f& left_p2 = left_line.end_point;
//		const Point2f& right_p1 = right_line.start_point;
//		const Point2f& right_p2 = right_line.end_point;
//
//		// ����ռ�ֱ�������˵�
//		const float f = 622.812; // ���࣬���㣬���߾���
//		const float cx = 654.193;
//		const float cy = 344.569;
//		const float Tx = 60.14;
//
//		float p1_disp = left_p1.x - right_p1.x;
//
//		float p1_z_w = f * Tx / p1_disp;
//		float p1_x_w = (left_p1.x - cx)*p1_z_w / f;
//		float p1_y_w = (left_p1.y - cy)*p1_z_w / f;
//
//		float p2_disp = left_p2.x - right_p2.x;
//		float p2_z_w = f * Tx / p2_disp;
//		float p2_x_w = (left_p2.x - cx)*p2_z_w / f;
//		float p2_y_w = (left_p2.y - cy)*p2_z_w / f;
//
//		// push�ؽ��Ŀռ�ֱ��
//		Line3D line_3d(Point3f(p1_x_w, p1_y_w, p1_z_w), Point3f(p2_x_w, p2_y_w, p2_z_w));
//		vec_line3d.push_back(line_3d);
//
//		// push �ռ�ֱ�߶˵�
//		pcl::PointXYZ p1_w(p1_x_w, p1_y_w, p1_z_w);
//		pcl::PointXYZ p2_w(p2_x_w, p2_y_w, p2_z_w);
//		cloud_reconstruct->points.push_back(p1_w);
//		cloud_reconstruct->points.push_back(p2_w);
//	}
//
//	// ȫ��ֱ�߲����ɵ���
//	float sample_dis = 1.0; // ֱ�ߴ�Ͼ���
//
//	for (const auto& line_3d : vec_line3d)
//	{
//		// ֱ�ߵȾ�����
//		float x0 = line_3d.x0_vec.x;
//		float y0 = line_3d.x0_vec.y;
//		float z0 = line_3d.x0_vec.z;
//
//		float v1 = line_3d.normal_vec.x;
//		float v2 = line_3d.normal_vec.y;
//		float v3 = line_3d.normal_vec.z;
//
//		float v_square = line_3d.normal_vec.dot(line_3d.normal_vec);
//		float deta_t = sample_dis / std::sqrt(v_square);
//
//		for (float t = line_3d.t_start; t < line_3d.t_end; t += deta_t)
//		{
//			pcl::PointXYZ point_3d;
//			point_3d.x = x0 + v1 * t;
//			point_3d.y = y0 + v2 * t;
//			point_3d.z = z0 + v3 * t;
//			cloud_reconstruct->points.push_back(point_3d);
//		}
//	}
//
//	//// ȫ��ֱ�߲����ɵ���
//	//float sample_dis = 1.0; // ֱ�ߴ�Ͼ���
//
//	//for (const auto& line_3d : vec_line3d)
//	//{
//	//	// ֱ�ߵȾ�����
//	//	Point3f deta_p = (line_3d.start_point - line_3d.end_point)/ line_3d.length;
//
//	//	for (float t = 0; t < line_3d.length; ++t)
//	//	{
//	//		static pcl::PointXYZ point_3d(line_3d.start_point.x, line_3d.start_point.y, line_3d.start_point.z);
//	//		point_3d.x = point_3d.x + deta_p.x;
//	//		point_3d.y = point_3d.y + deta_p.y;
//	//		point_3d.z = point_3d.z + deta_p.z;
//	//		cloud_reconstruct->points.push_back(point_3d);
//	//	}
//	//}
//
//	// ������׼
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src(new pcl::PointCloud<pcl::PointXYZ>());
//	createLWireFramePointCloud(cloud_src);
//	// ������Ʊ��淨��
//	clock_t fpfh_start_time = clock();
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_src; // Դ
//	ne_src.setInputCloud(cloud_src);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_src.setSearchMethod(tree_src);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_src.setRadiusSearch(2.0);
//	ne_src.setNumberOfThreads(8);
//	ne_src.compute(*cloud_src_normals);
//
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_tar; // �ؽ�
//	ne_tar.setInputCloud(cloud_reconstruct);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_reconstruct(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_tar.setSearchMethod(tree_reconstruct);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_reconstruct_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_tar.setRadiusSearch(2.0);
//	ne_src.setNumberOfThreads(8);
//	ne_tar.compute(*cloud_reconstruct_normals);
//
//	//����FPFH
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_src; // Դ
//	fpfh_src.setInputCloud(cloud_src);
//	fpfh_src.setInputNormals(cloud_src_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_src_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_src.setSearchMethod(tree_src_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_src(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_src.setRadiusSearch(1.8);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_src.compute(*fpfhs_src);
//	std::cout << "compute *cloud_src fpfh" << endl;
//
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_tar; // �ؽ�
//	fpfh_tar.setInputCloud(cloud_reconstruct);
//	fpfh_tar.setInputNormals(cloud_reconstruct_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_reconstruct_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_tar.setSearchMethod(tree_reconstruct_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_reconstruct(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_tar.setRadiusSearch(1.8);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_tar.compute(*fpfhs_reconstruct);
//
//	std::cout << "compute *cloud_tar fpfh" << endl;
//	clock_t fpfh_end_time = clock();
//	cout << "fpfh time: " << (double)(fpfh_end_time - fpfh_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//SAC��׼
//	clock_t sac_start_time = clock();
//
//	pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> scia;
//	scia.setInputSource(cloud_src);
//	scia.setInputTarget(cloud_reconstruct);
//	scia.setSourceFeatures(fpfhs_src);
//	scia.setTargetFeatures(fpfhs_reconstruct);
//
//	//scia.setMinSampleDistance(1);
//	//scia.setNumberOfSamples(2);
//	//scia.setCorrespondenceRandomness(20);
//
//	pcl::PointCloud<pcl::PointXYZ>::Ptr sac_result(new pcl::PointCloud<pcl::PointXYZ>);
//	scia.align(*sac_result);
//	Eigen::Matrix4f sac_trans = scia.getFinalTransformation();
//	clock_t sac_end_time = clock();
//	cout << "sac time: " << (double)(sac_end_time - sac_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//icp��׼
//	clock_t icp_start_time = clock();
//	pcl::PointCloud<pcl::PointXYZ>::Ptr icp_result(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
//	icp.setInputSource(cloud_src);
//	icp.setInputTarget(cloud_reconstruct);
//	// ��������Ӧ����� (e.g., correspondences with higher distances will be ignored)
//	icp.setMaxCorrespondenceDistance(40);
//	// ����������
//	icp.setMaximumIterations(100);
//	// ���α仯����֮��Ĳ�ֵ
//	icp.setTransformationEpsilon(1e-10);
//	// �������
//	icp.setEuclideanFitnessEpsilon(0.2);
//	icp.align(*icp_result, sac_trans);
//	clock_t icp_end_time = clock();
//	cout << "icp time: " << (double)(icp_end_time - icp_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//	// ��ȡ���յı任����
//	Eigen::Matrix4f final_trans_pcl = icp.getFinalTransformation();
//	cv::Mat_<float> final_trans_cv(4, 4);
//	for (size_t i = 0; i < 4; i++)
//	{
//		for (size_t j = 0; j < 4; j++)
//		{
//			final_trans_cv(i, j) = final_trans_pcl(i, j);
//		}
//	}
//
//	// ��������ϵ������
//	vector<cv::Point3f> vec_axis_3d;
//	vector<cv::Point2f> vec_axis_2d;
//	vec_axis_3d.push_back(Point3f(10.0, 0.0, 0.0));
//	vec_axis_3d.push_back(Point3f(0.0, 10.0, 0.0));
//	vec_axis_3d.push_back(Point3f(0.0, 0.0, 10.0));
//	vec_axis_3d.push_back(Point3f(0.0, 0.0, 0.0));
//
//	// ��������ϵת��������
//	Eigen::Vector3f final_angle;
//	Eigen::Vector3f final_translate;
//	matrix2Angle(final_trans_pcl, final_angle);
//	matrix2Translate(final_trans_pcl, final_translate);
//	cv::Mat rVec(3, 1, cv::DataType<float>::type); // Rotation vector
//	rVec.at<float>(0) = final_angle.x();
//	rVec.at<float>(1) = final_angle.y();
//	rVec.at<float>(2) = final_angle.z();
//
//	cv::Mat tVec(3, 1, cv::DataType<float>::type); // Translation vector
//	tVec.at<float>(0) = final_translate.x();
//	tVec.at<float>(1) = final_translate.y();
//	tVec.at<float>(2) = final_translate.z();
//
//	cv::projectPoints(vec_axis_3d, rVec, tVec, Pl.colRange(0, 3), Mat(), vec_axis_2d);
//
//	// ������ͶӰ���������
//	Mat reproject_axis_img = rleftimg.clone();
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[0], Scalar(0, 0, 255), 2);
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[1], Scalar(0, 255, 0), 2);
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[2], Scalar(255, 0, 0), 2);
//	imshow("reproject_axis_img", reproject_axis_img);
//	waitKey(0);
//
//	// ��ʾ����
//	// ��ɫ��Դ����
//	// ��ɫ��Ŀ�����
//	// ��ɫ����׼����
//	pcl::visualization::PCLVisualizer viewer("registration Viewer");
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h(cloud_src, 0, 255, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> tgt_h(cloud_reconstruct, 255, 0, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> final_h(icp_result, 0, 0, 255);
//	viewer.setBackgroundColor(1.0, 1.0, 1.0);
//
//	viewer.addPointCloud(cloud_src, src_h, "tar_points");
//	viewer.addPointCloud(cloud_reconstruct, tgt_h, "src_points");
//	viewer.addPointCloud(icp_result, final_h, "icp_points");
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "src_points"); // ���õ��ƴ�С
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "tar_points");
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "icp_points");
//	viewer.addCoordinateSystem(20.0); // �������ָʾ�ᣬ���ű���-------��ɫX����ɫY����ɫZ
//
//	//// �������
//	//cloud_reconstruct->width = 1;
//	//cloud_reconstruct->height = cloud_reconstruct->points.size();
//	//pcl::io::savePCDFileASCII("./cloud_reconstruct/cloud_reconstruct.pcd", *cloud_reconstruct);
//	//cout << "cloud_reconstruct is saved" << endl;
//
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(100);
//		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	}
//
//	return 0;
//}
//
//
//// L����
//void createLWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame)
//{
//
//	// 50
//	for (int i = 0; i <= 50; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i; p1.y = 0; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 0; p2.y = i; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//	}
//
//
//	// 35
//	for (int i = 0; i <= 35; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 14 + i; p1.y = 14; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 14; p2.y = 14 + i; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//	}
//
//	// 14
//	for (int i = 0; i <= 14; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 50; p1.y = i; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = i; p2.y = 50; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//	}
//}
//
//
////����תƽ�ƾ��������ת�Ƕ�
//void matrix2Angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle)
//{
//	double ax, ay, az;
//	if (result_trans(2, 0) == 1 || result_trans(2, 0) == -1)
//	{
//		az = 0;
//		double dlta;
//		dlta = atan2(result_trans(0, 1), result_trans(0, 2));
//		if (result_trans(2, 0) == -1)
//		{
//			ay = M_PI / 2;
//			ax = az + dlta;
//		}
//		else
//		{
//			ay = -M_PI / 2;
//			ax = -az + dlta;
//		}
//	}
//	else
//	{
//		ay = -asin(result_trans(2, 0));
//		ax = atan2(result_trans(2, 1) / cos(ay), result_trans(2, 2) / cos(ay));
//		az = atan2(result_trans(1, 0) / cos(ay), result_trans(0, 0) / cos(ay));
//	}
//	result_angle << ax, ay, az;
//}
//
////����תƽ�ƾ������ƽ�ƾ���
//void matrix2Translate(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_translate)
//{
//	result_translate << result_trans(0, 3), result_trans(1, 3), result_trans(2, 3);
//}




///*
//	L��
//*/
//
//#pragma once
//#include "LineAlgorithm.h"
//
//#include <iostream>  
//#include <time.h>
//
//#include <opencv2/opencv.hpp>
//#include <opencv2/ximgproc.hpp>
//#include <opencv2/xfeatures2d.hpp>
//
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/common/io.h>
//#include <pcl/point_cloud.h>
//
//#include <pcl/PolygonMesh.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <pcl/io/vtk_lib_io.h>
//#include <pcl/io/io.h>
//#include <pcl/io/pcd_io.h>//pcd ��д����ص�ͷ�ļ���
//#include <pcl/point_types.h> //PCL��֧�ֵĵ�����ͷ�ļ�
//#include <pcl/registration/ia_ransac.h>
//#include <pcl/point_cloud.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/features/normal_3d_omp.h>
//#include <pcl/features/fpfh.h>
//#include <pcl/features/fpfh_omp.h>
//#include <pcl/search/kdtree.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/filter.h>
//#include <pcl/registration/icp.h>
//#include <pcl/visualization/pcl_visualizer.h>
//
//using namespace cv;
//using namespace cv::ximgproc;
//using namespace std;
//
//// ͼƬ��С
//const int img_Width = 1280;
//const int img_Height = 720;
//
////У����ת����R��ͶӰ����P ��ͶӰ����Q
//Mat Rl, Rr, Pl, Pr, Q;
//
//// ����ģ�����
//void createLWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame);
//
//// ����ת��תƽ������
//void matrix2Angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle);
//void matrix2Translate(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_translate);
//
//int main()
//{
//
//	Mat left_and_right = imread("./left_right_img/IMG_0022.jpg");
//	Rect rect_left(0, 0, img_Width, img_Height);
//	Rect rect_right(img_Width, 0, img_Width, img_Height);
//
//	Mat left = left_and_right(rect_left);
//	Mat right = left_and_right(rect_right);
//
//	//imshow("src_left", left);
//	//imshow("src_right", right);
//	//waitKey(0);
//
//	Mat cameraMatrix[2], distCoeffs[2];
//	Mat Q;
//
//	// ����ڲ�
//	cameraMatrix[0] = Mat((Mat_<double>(3, 3) << 731.5022, 0.0, 0.0,
//		0.0, 732.0644, 0.0,
//		674.6773, 357.5894, 1.0)).t();
//
//	cameraMatrix[1] = Mat((Mat_<double>(3, 3) << 727.3993, 0.0, 0.0,
//		0.0, 727.6805, 0.0,
//		634.0939, 339.6110, 1.0)).t();
//
//	// ����ϵ��
//	distCoeffs[0] = (Mat_<double>(4, 1) << 0.1070, -0.1208, 0.0, -0.0004);
//	distCoeffs[1] = (Mat_<double>(4, 1) << 0.1145, -0.1517, 0.0007, -0.0004);
//
//	// ������λ��
//	Mat R = Mat((Mat_<double>(3, 3) << 1.0, 4.7372e-04, 4.8307e-04,
//		-4.7146e-04, 1.0, -0.0047,
//		-0.0005, 0.0047, 1.0)).t();
//	Mat T = (Mat_<double>(3, 1) << -60.1450, -0.0701, -0.0956);
//
//	// ͼ��У��֮�����Ч����
//	Rect validRoi[2];
//
//	// bougust����У��
//	stereoRectify(cameraMatrix[0], distCoeffs[0],
//		cameraMatrix[1], distCoeffs[1],
//		Size(img_Width, img_Height), R, T, Rl, Rr, Pl, Pr, Q,
//		CALIB_ZERO_DISPARITY, 1, Size(img_Width, img_Height), &validRoi[0], &validRoi[1]); //validPixROI1  validPixROI2 - ��ѡ�����������Rect�����ݡ����ڲ����������ض���Ч
//
//	//��ʾ���������ͶӰ����
//	cout << "camL" << cameraMatrix[0] << endl;
//	cout << "camR" << cameraMatrix[1] << endl;
//	cout << "Pro1" << Pl << endl;
//	cout << "Pro2" << Pr << endl;
//	cout << "Q" << Q << endl;
//
//	// ����������ͼ��У������ӳ���
//	Mat rmap[2][2];
//	initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], Rl, Pl, Size(img_Width, img_Height), CV_16SC2, rmap[0][0], rmap[0][1]);
//	initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], Rr, Pr, Size(img_Width, img_Height), CV_16SC2, rmap[1][0], rmap[1][1]);
//
//	// ��Դͼ�������ӳ�䣬ӳ�����ΪУ�����ͼ��
//	Mat rleftimg, rrightimg, rcolor;
//	remap(left, rleftimg, rmap[0][0], rmap[0][1], INTER_LINEAR);
//	remap(right, rrightimg, rmap[1][0], rmap[1][1], INTER_LINEAR);
//
//	////// ����У����ͼƬ
//	////imwrite("rleft.jpg", rleftimg);
//	////imwrite("rright.jpg", rrightimg);
//
//	//// �������߽��жԱ�
//	//Mat canvas;
//	//canvas.create(img_Height, img_Width * 2, CV_8UC3);
//	//hconcat(rleftimg, rrightimg, canvas);
//	//for (int j = 0; j < canvas.rows; j += 16)
//	//	line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
//	//line(canvas, Point(1280, 0), Point(1280, 720), Scalar(255, 255, 0), 1, 8);
//
//	//// ��ʾУ�����
//	//resize(canvas, canvas, Size(), 0.5, 0.5);
//	//namedWindow("rectifyResult", WINDOW_AUTOSIZE);
//	//imshow("rectifyResult", canvas);
//	//waitKey(0);
//
//
//	// �ҶȻ�
//	Mat left_img_gray, right_img_gray;
//	cvtColor(rleftimg, left_img_gray, COLOR_BGR2GRAY);
//	cvtColor(rrightimg, right_img_gray, COLOR_BGR2GRAY);
//	// ��
//	Mat left_img_sharpen, right_img_sharpen;
//	Mat kernel(3, 3, CV_32F, cv::Scalar(0));
//	kernel.at<float>(1, 1) = 5.0;
//	kernel.at<float>(0, 1) = -1.0;
//	kernel.at<float>(1, 0) = -1.0;
//	kernel.at<float>(1, 2) = -1.0;
//	kernel.at<float>(2, 1) = -1.0;
//	filter2D(left_img_gray, left_img_sharpen, left_img_gray.depth(), kernel);
//	filter2D(right_img_gray, right_img_sharpen, right_img_gray.depth(), kernel);
//
//	// YOLOv3 Ŀ����
//	Rect rect_left_tar(850, 200, 240, 240);
//	Rect rect_right_tar(650, 200, 270, 240);
//	Mat left_tar_img = left_img_sharpen(rect_left_tar); // ��ȡĿ������
//	Mat right_tar_img = right_img_sharpen(rect_right_tar);
//
//	 //����Ŀ������
//	//imshow("left_lines", left_tar_img);
//	//imshow("right_lines", right_tar_img);
//	///imwrite("roi1.jpg", left_tar_img);
//	///imwrite("roi2.jpg", right_tar_img);
//	//waitKey(0);
//
//	// ��Ŀ���������ֱ�߼��
//	int    length_threshold = 20; // ���ڴ���ֵ��ֱ�߱�����
//	float  distance_threshold = 1.41421356f;  // ������߶θ�Զ�ĵ㱻��Ϊ����Ⱥ��
//	double canny_th1 = 40.0; // �ݶ�С�ڴ���ֵ����Ϊ�Ǳ�Ե
//	double canny_th2 = 100.0; // �ݶȴ��ڴ���ֵһ���Ǳ�Ե
//	int    canny_aperture_size = 3;
//	bool   do_merge = true;
//	Ptr<FastLineDetector> fld = createFastLineDetector(
//		length_threshold,
//		distance_threshold,
//		canny_th1,
//		canny_th2,
//		canny_aperture_size,
//		do_merge);
//
//	vector<Vec4f> left_lines_fld;
//	vector<Vec4f> right_lines_fld;
//	fld->detect(left_tar_img, left_lines_fld);
//	fld->detect(right_tar_img, right_lines_fld);
//
//	//// ����fld�����
//	//fld->drawSegments(left_tar_img, left_lines_fld, false);
//	//fld->drawSegments(right_tar_img, right_lines_fld, false);
//	//imshow("left_lines", left_tar_img);
//	//imshow("right_lines", right_tar_img);
//	//waitKey(0);
//
//	// ֱ�߸�ʽת��
//	vector<Line2D> vec_left_tar_line2d;
//	vector<Line2D> vec_right_tar_line2d;
//	for (const auto& line_4f : left_lines_fld) //��
//	{
//		static int serial_number = 0;
//		if (serial_number == 0)
//		{
//			vec_left_tar_line2d.push_back(Line2D(serial_number, Vec4f(164.0,107.0,85.84,10.32)));
//			++serial_number;
//			continue;
//		}
//		if (serial_number == 1)
//		{
//			vec_left_tar_line2d.push_back(Line2D(serial_number, Vec4f(74, 14, 38.019, 116.689)));
//			++serial_number;
//			continue;
//		}
//		if (serial_number == 4)
//		{
//			++serial_number;
//			continue;
//		}
//		if (serial_number == 6)
//		{
//			vec_left_tar_line2d.push_back(Line2D(serial_number, Vec4f(99, 212, 149.9433, 110.9703)));
//			++serial_number;
//			continue;
//		}
//		if (serial_number == 9)
//		{
//			vec_left_tar_line2d.push_back(Line2D(serial_number, Vec4f(100, 212, 224, 192)));
//			++serial_number;
//			continue;
//		}
//		vec_left_tar_line2d.push_back(Line2D(serial_number, line_4f));
//		++serial_number;
//	}
//	vec_left_tar_line2d.insert(vec_left_tar_line2d.begin() + 4, Line2D(4, Vec4f(75, 15, 86, 11)));   //��ָ��λ�ã������ڵ����Ԫ��ǰ����һ��Ԫ��
//	vec_left_tar_line2d.push_back(Line2D(10, Vec4f(151, 111, 164, 109)));
//	vec_left_tar_line2d.push_back(Line2D(11, Vec4f(230, 183, 225, 192)));
//
//
//	for (const auto& line_4f : right_lines_fld) // ��
//	{
//		static int serial_number = 0;
//		if (serial_number == 7)
//		{
//			vec_right_tar_line2d.push_back(Line2D(serial_number, Vec4f(84, 210, 55, 118)));
//			++serial_number;
//			continue;
//		}
//		if (serial_number == 8)
//		{
//			vec_right_tar_line2d.push_back(Line2D(serial_number, Vec4f(103.114, 200.729, 223, 183)));
//			++serial_number;
//			continue;
//		}
//		if (serial_number == 9)
//		{
//			vec_right_tar_line2d.push_back(Line2D(serial_number, Vec4f(88, 211, 221, 191)));
//			++serial_number;
//			continue;
//		}
//
//		vec_right_tar_line2d.push_back(Line2D(serial_number, line_4f));
//		++serial_number;
//	}
//	
//	vec_right_tar_line2d.push_back(Line2D(10, Vec4f(98, 110, 112, 109)));
//	vec_right_tar_line2d.push_back(Line2D(11, Vec4f(60, 13, 72, 12)));
//	vec_right_tar_line2d.push_back(Line2D(12, Vec4f(221, 191, 222, 183)));
//
//	//// ����ֱ�ߺ�������
//	//Mat left_line_img = left_tar_img.clone();
//	//Mat right_line_img = right_tar_img.clone();
//	//cvtColor(left_line_img, left_line_img, COLOR_GRAY2BGR);
//	//cvtColor(right_line_img, right_line_img, COLOR_GRAY2BGR);
//	//for (const auto& line2d : vec_left_tar_line2d)
//	//{
//	//	cv::line(left_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(left_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_tar_line2d)
//	//{
//	//	cv::line(right_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(right_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//imshow("left_lines", left_line_img);
//	//imshow("right_lines", right_line_img);
//	//waitKey(0);
//
//	// ��ֱ�߶˵�ȡ����
//	vector<int> left_line_points_index;
//	vector<int> right_line_points_index;
//	vector<LinePoint2d> vec_left_line_points; // ������ֱ�߶˵�
//	vector<LinePoint2d> vec_right_line_points; // ������ֱ�߶˵�
//	for (const auto& line2d : vec_left_tar_line2d) // ��
//	{
//		static int index = 0;
//
//		LinePoint2d line_point1(line2d.start_point, index, line2d.serial_number, true);
//		LinePoint2d line_point2(line2d.end_point, index + 1, line2d.serial_number, false);
//
//		vec_left_line_points.push_back(line_point1);
//		vec_left_line_points.push_back(line_point2);
//		left_line_points_index.push_back(index);
//		left_line_points_index.push_back(index + 1);
//
//		index += 2;
//	}
//	for (const auto& line2d : vec_right_tar_line2d) // ��
//	{
//		static int index = 0;
//
//		LinePoint2d line_point1(line2d.start_point, index, line2d.serial_number, true);
//		LinePoint2d line_point2(line2d.end_point, index + 1, line2d.serial_number, false);
//
//		vec_right_line_points.push_back(line_point1);
//		vec_right_line_points.push_back(line_point2);
//		right_line_points_index.push_back(index);
//		right_line_points_index.push_back(index + 1);
//
//		index += 2;
//	}
//
//	// �Զ˵������ֱ�߽��з���
//	vector<vector<int>> vec_left_near_points_group; // �洢����vec_left_line_points���±�
//	vector<vector<int>> vec_right_near_points_group;
//
//	double line_point_near_threshold = 10; // �˵��������ֵ
//
//	while (left_line_points_index.size() != 0)
//	{
//		int curr_index = left_line_points_index[0];
//		vector<int> curr_near;
//		LinePoint2d curr_point = vec_left_line_points[curr_index];
//		for (int j = 0; j < left_line_points_index.size(); j++)
//		{
//			int comp_index = left_line_points_index[j];
//
//			LinePoint2d comp_point = vec_left_line_points[comp_index];
//
//			double dis = computeTowPointDistance(curr_point.point_2d, comp_point.point_2d);
//
//			if (dis <= line_point_near_threshold)
//			{
//				curr_near.push_back(comp_index);
//			}
//		}
//
//		// ɾ���Ѿ���ԵĶ˵�
//		for (int i = 0; i < curr_near.size(); i++)
//		{
//			deleteElement(left_line_points_index, curr_near[i]);
//		}
//		vec_left_near_points_group.push_back(curr_near);
//	}
//
//	while (right_line_points_index.size() != 0)
//	{
//		int curr_index = right_line_points_index[0];
//		vector<int> curr_near;
//		LinePoint2d curr_point = vec_right_line_points[curr_index];
//		for (int j = 0; j < right_line_points_index.size(); j++)
//		{
//			int comp_index = right_line_points_index[j];
//
//			LinePoint2d comp_point = vec_right_line_points[comp_index];
//
//			double dis = computeTowPointDistance(curr_point.point_2d, comp_point.point_2d);
//
//			if (dis <= line_point_near_threshold)
//			{
//				curr_near.push_back(comp_index);
//			}
//		}
//
//		// ɾ���Ѿ���ԵĶ˵�
//		for (int i = 0; i < curr_near.size(); i++)
//		{
//			deleteElement(right_line_points_index, curr_near[i]);
//		}
//		vec_right_near_points_group.push_back(curr_near);
//	}
//
//	// ֱ�ߺϲ�
//	double line_slope_similarity_threshold = 0.5; // ֱ��б�����ƶ���ֵ
//
//	// �Զ˵������ֱ�߽����ӳ����󽻵�
//	Mat left_cross_point_image = left_tar_img.clone();
//	Mat right_cross_point_image = right_tar_img.clone();
//	cvtColor(left_cross_point_image, left_cross_point_image, COLOR_GRAY2BGR);
//	cvtColor(right_cross_point_image, right_cross_point_image, COLOR_GRAY2BGR);
//	for (const auto& near_points : vec_left_near_points_group) // ��
//	{
//		if (near_points.size() == 1)
//			continue;
//		if (near_points.size() == 2)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			LinePoint2d point1 = vec_left_line_points[point_index1];
//			LinePoint2d point2 = vec_left_line_points[point_index2];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			Line2D& line1 = vec_left_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_left_tar_line2d[line_serial_number2];
//			Point2f cross_point = conputeIntersectionPoint(line1, line2);
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			cv::circle(left_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//		if (near_points.size() == 3)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			int point_index3 = near_points[2];
//			LinePoint2d point1 = vec_left_line_points[point_index1];
//			LinePoint2d point2 = vec_left_line_points[point_index2];
//			LinePoint2d point3 = vec_left_line_points[point_index3];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			int line_serial_number3 = point3.line_serial_number;
//			Line2D& line1 = vec_left_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_left_tar_line2d[line_serial_number2];
//			Line2D& line3 = vec_left_tar_line2d[line_serial_number3];
//			Point2f cross_point1 = conputeIntersectionPoint(line1, line2);
//			Point2f cross_point2 = conputeIntersectionPoint(line1, line3);
//			Point2f cross_point3 = conputeIntersectionPoint(line2, line3);
//			Point2f cross_point = (cross_point1 + cross_point2 + cross_point3) / 3.0;
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			if (point3.is_start_point == true)
//				line3.start_point = cross_point;
//			else
//				line3.end_point = cross_point;
//
//			cv::circle(left_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//	}
//	for (const auto& near_points : vec_right_near_points_group) // ��
//	{
//		if (near_points.size() == 1)
//			continue;
//		if (near_points.size() == 2)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			LinePoint2d point1 = vec_right_line_points[point_index1];
//			LinePoint2d point2 = vec_right_line_points[point_index2];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			Line2D& line1 = vec_right_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_right_tar_line2d[line_serial_number2];
//			Point2f cross_point = conputeIntersectionPoint(line1, line2);
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			cv::circle(right_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//		if (near_points.size() == 3)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			int point_index3 = near_points[2];
//			LinePoint2d point1 = vec_right_line_points[point_index1];
//			LinePoint2d point2 = vec_right_line_points[point_index2];
//			LinePoint2d point3 = vec_right_line_points[point_index3];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			int line_serial_number3 = point3.line_serial_number;
//			Line2D& line1 = vec_right_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_right_tar_line2d[line_serial_number2];
//			Line2D& line3 = vec_right_tar_line2d[line_serial_number3];
//			Point2f cross_point1 = conputeIntersectionPoint(line1, line2);
//			Point2f cross_point2 = conputeIntersectionPoint(line1, line3);
//			Point2f cross_point3 = conputeIntersectionPoint(line2, line3);
//			Point2f cross_point = (cross_point1 + cross_point2 + cross_point3) / 3.0;
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			if (point3.is_start_point == true)
//				line3.start_point = cross_point;
//			else
//				line3.end_point = cross_point;
//
//			cv::circle(right_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//	}
//
//	//// ��ʾֱ�߽���
//	//imshow("left_cross_point_image", left_cross_point_image);
//	//imshow("right_cross_point_image", right_cross_point_image);
//	//waitKey(0);
//
//	//// ��ʾ����������Ϊ�˵��ֱ��
//	//Mat left_line_img = left_tar_img.clone();
//	//Mat right_line_img = right_tar_img.clone();
//	//cvtColor(left_line_img, left_line_img, COLOR_GRAY2BGR);
//	//cvtColor(right_line_img, right_line_img, COLOR_GRAY2BGR);
//	//for (const auto& line2d : vec_left_tar_line2d)
//	//{
//	//	cv::line(left_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(left_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_tar_line2d)
//	//{
//	//	cv::line(right_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(right_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//imshow("left_lines", left_line_img);
//	//imshow("right_lines", right_line_img);
//	//waitKey(0);
//
//
//	// ֱ�߾ֲ�����ӳ�䵽ȫ�����꣬���¼���б��
//	vector<Line2D> vec_left_line2d;
//	vector<Line2D> vec_right_line2d;
//	for (const auto& line2d : vec_left_tar_line2d) // ��
//	{
//		static int index = 0;
//		float x1 = line2d.start_point.x + rect_left_tar.x;
//		float y1 = line2d.start_point.y + rect_left_tar.y;
//		float x2 = line2d.end_point.x + rect_left_tar.x;
//		float y2 = line2d.end_point.y + rect_left_tar.y;
//		Vec4f line4f(x1, y1, x2, y2);
//		vec_left_line2d.push_back(Line2D(index, line4f));
//		++index;
//	}
//	for (const auto& line2d : vec_right_tar_line2d) // ��
//	{
//		static int index = 0;
//		float x1 = line2d.start_point.x + rect_right_tar.x;
//		float y1 = line2d.start_point.y + rect_right_tar.y;
//		float x2 = line2d.end_point.x + rect_right_tar.x;
//		float y2 = line2d.end_point.y + rect_right_tar.y;
//		Vec4f line4f(x1, y1, x2, y2);
//		vec_right_line2d.push_back(Line2D(index, line4f));
//		++index;
//	}
//
//	//// ���ƴ���õ�ȫ��ֱ��
//	//Mat left_global_line_img = rleftimg.clone();
//	//Mat right_global_line_img = rrightimg.clone();
//
//	//for (const auto& line2d : vec_left_line2d)
//	//{
//	//	cv::line(left_global_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 2);
//	//	cv::putText(left_global_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_line2d)
//	//{
//	//	cv::line(right_global_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 2);
//	//	cv::putText(right_global_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//
//	//Mat left_right_global_line_img; // ƴ������ͼ��
//	//left_right_global_line_img.create(img_Height, img_Width * 2, CV_8UC3);
//	//hconcat(left_global_line_img, right_global_line_img, left_right_global_line_img);
//	//imshow("left_global_line_img", left_global_line_img);
//	//imshow("right_global_line_img", right_global_line_img);
//	//imshow("left_right_global_line_img", left_right_global_line_img);
//	//waitKey(0);
//
//	// ֱ��ƥ��
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_reconstruct(new pcl::PointCloud<pcl::PointXYZ>()); // ֱ�ߵ���
//	std::vector<std::vector<int>> vec_match_index_group; // ƥ������
//	//lineMatch(vec_left_line2d, vec_right_line2d, vec_match_index_group);
//	std::vector<int> match1; match1.push_back(0); match1.push_back(0); vec_match_index_group.push_back(match1);
//	std::vector<int> match2; match2.push_back(1); match2.push_back(2); vec_match_index_group.push_back(match2);
//	std::vector<int> match3; match3.push_back(2); match3.push_back(1); vec_match_index_group.push_back(match3);
//	std::vector<int> match4; match4.push_back(3); match4.push_back(4); vec_match_index_group.push_back(match4);
//	//std::vector<int> match5; match5.push_back(4); match5.push_back(11); vec_match_index_group.push_back(match5);
//	std::vector<int> match6; match6.push_back(5); match6.push_back(5); vec_match_index_group.push_back(match6);
//	std::vector<int> match7; match7.push_back(6); match7.push_back(6); vec_match_index_group.push_back(match7);
//	std::vector<int> match8; match8.push_back(7); match8.push_back(7); vec_match_index_group.push_back(match8);
//	std::vector<int> match9; match9.push_back(8); match9.push_back(8); vec_match_index_group.push_back(match9);
//	std::vector<int> match10; match10.push_back(9); match10.push_back(9); vec_match_index_group.push_back(match10);
//	std::vector<int> match11; match11.push_back(10); match11.push_back(10); vec_match_index_group.push_back(match11);
//	std::vector<int> match12; match12.push_back(11); match12.push_back(12); vec_match_index_group.push_back(match12);
//
//	// ֱ���ؽ�
//	vector<Line3D> vec_line3d;
//	for (const auto vec_match_index : vec_match_index_group)
//	{
//		if (vec_match_index.size() != 2)
//			continue;
//
//		int left_index = vec_match_index[0];
//		int right_index = vec_match_index[1];
//		const Line2D& left_line = vec_left_line2d[left_index];
//		const Line2D& right_line = vec_right_line2d[right_index];
//
//		// p1��㣬p2�յ�
//		const Point2f& left_p1 = left_line.start_point;
//		const Point2f& left_p2 = left_line.end_point;
//		const Point2f& right_p1 = right_line.start_point;
//		const Point2f& right_p2 = right_line.end_point;
//
//		// ����ռ�ֱ�������˵�
//		const float f = 622.812; // ���࣬���㣬���߾���
//		const float cx = 654.193;
//		const float cy = 344.569;
//		const float Tx = 60.14;
//
//		float p1_disp = left_p1.x - right_p1.x;
//
//		float p1_z_w = f * Tx / p1_disp;
//		float p1_x_w = (left_p1.x - cx)*p1_z_w / f;
//		float p1_y_w = (left_p1.y - cy)*p1_z_w / f;
//
//		float p2_disp = left_p2.x - right_p2.x;
//		float p2_z_w = f * Tx / p2_disp;
//		float p2_x_w = (left_p2.x - cx)*p2_z_w / f;
//		float p2_y_w = (left_p2.y - cy)*p2_z_w / f;
//
//		// push�ؽ��Ŀռ�ֱ��
//		Line3D line_3d(Point3f(p1_x_w, p1_y_w, p1_z_w), Point3f(p2_x_w, p2_y_w, p2_z_w));
//		vec_line3d.push_back(line_3d);
//
//		// push �ռ�ֱ�߶˵�
//		pcl::PointXYZ p1_w(p1_x_w, p1_y_w, p1_z_w);
//		pcl::PointXYZ p2_w(p2_x_w, p2_y_w, p2_z_w);
//		cloud_reconstruct->points.push_back(p1_w);
//		cloud_reconstruct->points.push_back(p2_w);
//	}
//
//	// ȫ��ֱ�߲����ɵ���
//	float sample_dis = 1.0; // ֱ�ߴ�Ͼ���
//
//	for (const auto& line_3d : vec_line3d)
//	{
//		// ֱ�ߵȾ�����
//		float x0 = line_3d.x0_vec.x;
//		float y0 = line_3d.x0_vec.y;
//		float z0 = line_3d.x0_vec.z;
//
//		float v1 = line_3d.normal_vec.x;
//		float v2 = line_3d.normal_vec.y;
//		float v3 = line_3d.normal_vec.z;
//
//		float v_square = line_3d.normal_vec.dot(line_3d.normal_vec);
//		float deta_t = sample_dis / std::sqrt(v_square);
//
//		for (float t = line_3d.t_start; t < line_3d.t_end; t += deta_t)
//		{
//			pcl::PointXYZ point_3d;
//			point_3d.x = x0 + v1 * t;
//			point_3d.y = y0 + v2 * t;
//			point_3d.z = z0 + v3 * t;
//			cloud_reconstruct->points.push_back(point_3d);
//		}
//	}
//
//	// ������׼
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src(new pcl::PointCloud<pcl::PointXYZ>());
//	createLWireFramePointCloud(cloud_src);
//	// ������Ʊ��淨��
//	clock_t fpfh_start_time = clock();
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_src; // Դ
//	ne_src.setInputCloud(cloud_src);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_src.setSearchMethod(tree_src);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_src.setRadiusSearch(2.0);
//	ne_src.setNumberOfThreads(8);
//	ne_src.compute(*cloud_src_normals);
//
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_tar; // �ؽ�
//	ne_tar.setInputCloud(cloud_reconstruct);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_reconstruct(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_tar.setSearchMethod(tree_reconstruct);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_reconstruct_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_tar.setRadiusSearch(2.0);
//	ne_src.setNumberOfThreads(8);
//	ne_tar.compute(*cloud_reconstruct_normals);
//
//	//����FPFH
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_src; // Դ
//	fpfh_src.setInputCloud(cloud_src);
//	fpfh_src.setInputNormals(cloud_src_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_src_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_src.setSearchMethod(tree_src_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_src(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_src.setRadiusSearch(1.8);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_src.compute(*fpfhs_src);
//	std::cout << "compute *cloud_src fpfh" << endl;
//
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_tar; // �ؽ�
//	fpfh_tar.setInputCloud(cloud_reconstruct);
//	fpfh_tar.setInputNormals(cloud_reconstruct_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_reconstruct_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_tar.setSearchMethod(tree_reconstruct_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_reconstruct(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_tar.setRadiusSearch(1.8);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_tar.compute(*fpfhs_reconstruct);
//
//	std::cout << "compute *cloud_tar fpfh" << endl;
//	clock_t fpfh_end_time = clock();
//	cout << "fpfh time: " << (double)(fpfh_end_time - fpfh_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//SAC��׼
//	clock_t sac_start_time = clock();
//
//	pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> scia;
//	scia.setInputSource(cloud_src);
//	scia.setInputTarget(cloud_reconstruct);
//	scia.setSourceFeatures(fpfhs_src);
//	scia.setTargetFeatures(fpfhs_reconstruct);
//
//	//scia.setMinSampleDistance(1);
//	//scia.setNumberOfSamples(2);
//	//scia.setCorrespondenceRandomness(20);
//
//	pcl::PointCloud<pcl::PointXYZ>::Ptr sac_result(new pcl::PointCloud<pcl::PointXYZ>);
//	scia.align(*sac_result);
//	Eigen::Matrix4f sac_trans = scia.getFinalTransformation();
//	clock_t sac_end_time = clock();
//	cout << "sac time: " << (double)(sac_end_time - sac_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//icp��׼
//	clock_t icp_start_time = clock();
//	pcl::PointCloud<pcl::PointXYZ>::Ptr icp_result(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
//	icp.setInputSource(cloud_src);
//	icp.setInputTarget(cloud_reconstruct);
//	// ��������Ӧ����� (e.g., correspondences with higher distances will be ignored)
//	icp.setMaxCorrespondenceDistance(40);
//	// ����������
//	icp.setMaximumIterations(100);
//	// ���α仯����֮��Ĳ�ֵ
//	icp.setTransformationEpsilon(1e-10);
//	// �������
//	icp.setEuclideanFitnessEpsilon(0.2);
//	icp.align(*icp_result, sac_trans);
//	clock_t icp_end_time = clock();
//	cout << "icp time: " << (double)(icp_end_time - icp_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//	// ��ȡ���յı任����
//	Eigen::Matrix4f final_trans_pcl = icp.getFinalTransformation();
//	cv::Mat_<float> final_trans_cv(4, 4);
//	for (size_t i = 0; i < 4; i++)
//	{
//		for (size_t j = 0; j < 4; j++)
//		{
//			final_trans_cv(i, j) = final_trans_pcl(i, j);
//		}
//	}
//
//	// ��������ϵ������
//	vector<cv::Point3f> vec_axis_3d;
//	vector<cv::Point2f> vec_axis_2d;
//	vec_axis_3d.push_back(Point3f(10.0, 0.0, 0.0));
//	vec_axis_3d.push_back(Point3f(0.0, 10.0, 0.0));
//	vec_axis_3d.push_back(Point3f(0.0, 0.0, 10.0));
//	vec_axis_3d.push_back(Point3f(0.0, 0.0, 0.0));
//
//	// ��������ϵת��������
//	Eigen::Vector3f final_angle;
//	Eigen::Vector3f final_translate;
//	matrix2Angle(final_trans_pcl, final_angle);
//	matrix2Translate(final_trans_pcl, final_translate);
//	cv::Mat rVec(3, 1, cv::DataType<float>::type); // Rotation vector
//	rVec.at<float>(0) = final_angle.x();
//	rVec.at<float>(1) = final_angle.y();
//	rVec.at<float>(2) = final_angle.z();
//
//	cv::Mat tVec(3, 1, cv::DataType<float>::type); // Translation vector
//	tVec.at<float>(0) = final_translate.x();
//	tVec.at<float>(1) = final_translate.y();
//	tVec.at<float>(2) = final_translate.z();
//
//	cout << "λ�ƣ�" << tVec << endl;
//	cout << "��̬��" << rVec << endl;
//
//	cv::projectPoints(vec_axis_3d, rVec, tVec, Pl.colRange(0, 3), Mat(), vec_axis_2d);
//
//	// ������ͶӰ���������
//	Mat reproject_axis_img = rleftimg.clone();
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[0], Scalar(0, 0, 255), 2);
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[1], Scalar(0, 255, 0), 2);
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[2], Scalar(255, 0, 0), 2);
//	imshow("reproject_axis_img", reproject_axis_img);
//	waitKey(0);
//
//	// ��ʾ����
//	// ��ɫ��Դ����
//	// ��ɫ��Ŀ�����
//	// ��ɫ����׼����
//	pcl::visualization::PCLVisualizer viewer("registration Viewer");
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h(cloud_src, 0, 255, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> tgt_h(cloud_reconstruct, 255, 0, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> final_h(icp_result, 0, 0, 255);
//	viewer.setBackgroundColor(1.0, 1.0, 1.0);
//
//	viewer.addPointCloud(cloud_src, src_h, "tar_points");
//	viewer.addPointCloud(cloud_reconstruct, tgt_h, "src_points");
//	viewer.addPointCloud(icp_result, final_h, "icp_points");
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "src_points"); // ���õ��ƴ�С
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "tar_points");
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "icp_points");
//	viewer.addCoordinateSystem(20.0); // �������ָʾ�ᣬ���ű���-------��ɫX����ɫY����ɫZ
//
//	// �����ؽ�����
//	cloud_reconstruct->width = 1;
//	cloud_reconstruct->height = cloud_reconstruct->points.size();
//	pcl::io::savePCDFileASCII("./cloud_reconstruct_src/l_reconstruct.pcd", *cloud_reconstruct);
//	cout << "cloud_reconstruct is saved" << endl;
//
//	// ����Դ����
//	cloud_src->width = 1;
//	cloud_src->height = cloud_src->points.size();
//	pcl::io::savePCDFileASCII("./cloud_reconstruct_src/l_src.pcd", *cloud_src);
//	cout << "cloud_src is saved" << endl;
//
//	// ������׼����
//	icp_result->width = 1;
//	icp_result->height = icp_result->points.size();
//	pcl::io::savePCDFileASCII("./cloud_reconstruct_src/l_icp_result.pcd", *icp_result);
//	cout << "icp_result is saved" << endl;
//
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(100);
//		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	}
//
//	return 0;
//}
//
//
//// L����
//void createLWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame)
//{
//	// 40
//	for (int i = 0; i <= 40; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i; p1.y = 0; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 0; p2.y = i; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 0; p3.y = 0; p3.z = i;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 40; p4.y = i; p4.z = 0;
//		cloud_wire_frame->points.push_back(p4);
//
//		pcl::PointXYZ p5;
//		p5.x = i; p5.y = 40; p5.z = 0;
//		cloud_wire_frame->points.push_back(p5);
//
//		pcl::PointXYZ p6;
//		p6.x = 0; p6.y = 40; p6.z = i;
//		cloud_wire_frame->points.push_back(p6);
//
//		pcl::PointXYZ p7;
//		p7.x = 0; p7.y = i; p7.z = 40;
//		cloud_wire_frame->points.push_back(p7);
//
//		pcl::PointXYZ p8;
//		p8.x = 5; p8.y = i; p8.z = 40;
//		cloud_wire_frame->points.push_back(p8);
//
//		pcl::PointXYZ p9;
//		p9.x = 40; p9.y = i; p9.z = 5;
//		cloud_wire_frame->points.push_back(p9);
//	}
//
//	// 5
//	for (int i = 0; i <= 5; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i; p1.y = 0; p1.z = 40;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = i; p2.y = 40; p2.z = 40;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 40; p3.y = 0; p3.z = i;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 40; p4.y = 40; p4.z = i;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//	// 35
//	for (int i = 0; i <= 35; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 5 + i; p1.y = 0; p1.z = 5;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 5 + i; p2.y = 40; p2.z = 5;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 5; p3.y = 0; p3.z = i + 5;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 5; p4.y = 40; p4.z = i + 5;
//		cloud_wire_frame->points.push_back(p4);
//	}
//	
//}
//
//
////����תƽ�ƾ��������ת�Ƕ�
//void matrix2Angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle)
//{
//	double ax, ay, az;
//	if (result_trans(2, 0) == 1 || result_trans(2, 0) == -1)
//	{
//		az = 0;
//		double dlta;
//		dlta = atan2(result_trans(0, 1), result_trans(0, 2));
//		if (result_trans(2, 0) == -1)
//		{
//			ay = M_PI / 2;
//			ax = az + dlta;
//		}
//		else
//		{
//			ay = -M_PI / 2;
//			ax = -az + dlta;
//		}
//	}
//	else
//	{
//		ay = -asin(result_trans(2, 0));
//		ax = atan2(result_trans(2, 1) / cos(ay), result_trans(2, 2) / cos(ay));
//		az = atan2(result_trans(1, 0) / cos(ay), result_trans(0, 0) / cos(ay));
//	}
//	result_angle << ax, ay, az;
//}
//
////����תƽ�ƾ������ƽ�ƾ���
//void matrix2Translate(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_translate)
//{
//	result_translate << result_trans(0, 3), result_trans(1, 3), result_trans(2, 3);
//}


///*
//  T��Ƭ
//*/
//
//#pragma once
//#include "LineAlgorithm.h"
//
//#include <iostream>  
//#include <time.h>
//
//#include <opencv2/opencv.hpp>
//#include <opencv2/ximgproc.hpp>
//#include <opencv2/xfeatures2d.hpp>
//
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/common/io.h>
//#include <pcl/point_cloud.h>
//
//#include <pcl/PolygonMesh.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <pcl/io/vtk_lib_io.h>
//#include <pcl/io/io.h>
//#include <pcl/io/pcd_io.h>//pcd ��д����ص�ͷ�ļ���
//#include <pcl/point_types.h> //PCL��֧�ֵĵ�����ͷ�ļ�
//#include <pcl/registration/ia_ransac.h>
//#include <pcl/point_cloud.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/features/normal_3d_omp.h>
//#include <pcl/features/fpfh.h>
//#include <pcl/features/fpfh_omp.h>
//#include <pcl/search/kdtree.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/filter.h>
//#include <pcl/registration/icp.h>
//#include <pcl/visualization/pcl_visualizer.h>
//
//using namespace cv;
//using namespace cv::ximgproc;
//using namespace std;
//
//// ͼƬ��С
//const int img_Width = 1280;
//const int img_Height = 720;
//
////У����ת����R��ͶӰ����P ��ͶӰ����Q
//Mat Rl, Rr, Pl, Pr, Q;
//
//// ����ģ�����
//void createTWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame);
//
//// ����ת��תƽ������
//void matrix2Angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle);
//void matrix2Translate(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_translate);
//
//int main()
//{
//
//	Mat left_and_right = imread("./left_right_img/IMG_0022.jpg");
//	Rect rect_left(0, 0, img_Width, img_Height);
//	Rect rect_right(img_Width, 0, img_Width, img_Height);
//
//	Mat left = left_and_right(rect_left);
//	Mat right = left_and_right(rect_right);
//
//	//imshow("src_left", left);
//	//imshow("src_right", right);
//	//waitKey(0);
//
//	Mat cameraMatrix[2], distCoeffs[2];
//	Mat Q;
//
//	// ����ڲ�
//	cameraMatrix[0] = Mat((Mat_<double>(3, 3) << 731.5022, 0.0, 0.0,
//		0.0, 732.0644, 0.0,
//		674.6773, 357.5894, 1.0)).t();
//
//	cameraMatrix[1] = Mat((Mat_<double>(3, 3) << 727.3993, 0.0, 0.0,
//		0.0, 727.6805, 0.0,
//		634.0939, 339.6110, 1.0)).t();
//
//	// ����ϵ��
//	distCoeffs[0] = (Mat_<double>(4, 1) << 0.1070, -0.1208, 0.0, -0.0004);
//	distCoeffs[1] = (Mat_<double>(4, 1) << 0.1145, -0.1517, 0.0007, -0.0004);
//
//	// ������λ��
//	Mat R = Mat((Mat_<double>(3, 3) << 1.0, 4.7372e-04, 4.8307e-04,
//		-4.7146e-04, 1.0, -0.0047,
//		-0.0005, 0.0047, 1.0)).t();
//	Mat T = (Mat_<double>(3, 1) << -60.1450, -0.0701, -0.0956);
//
//	// ͼ��У��֮�����Ч����
//	Rect validRoi[2];
//
//	// bougust����У��
//	stereoRectify(cameraMatrix[0], distCoeffs[0],
//		cameraMatrix[1], distCoeffs[1],
//		Size(img_Width, img_Height), R, T, Rl, Rr, Pl, Pr, Q,
//		CALIB_ZERO_DISPARITY, 1, Size(img_Width, img_Height), &validRoi[0], &validRoi[1]); //validPixROI1  validPixROI2 - ��ѡ�����������Rect�����ݡ����ڲ����������ض���Ч
//
//	//��ʾ���������ͶӰ����
//	cout << "camL" << cameraMatrix[0] << endl;
//	cout << "camR" << cameraMatrix[1] << endl;
//	cout << "Pro1" << Pl << endl;
//	cout << "Pro2" << Pr << endl;
//	cout << "Q" << Q << endl;
//
//	// ����������ͼ��У������ӳ���
//	Mat rmap[2][2];
//	initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], Rl, Pl, Size(img_Width, img_Height), CV_16SC2, rmap[0][0], rmap[0][1]);
//	initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], Rr, Pr, Size(img_Width, img_Height), CV_16SC2, rmap[1][0], rmap[1][1]);
//
//	// ��Դͼ�������ӳ�䣬ӳ�����ΪУ�����ͼ��
//	Mat rleftimg, rrightimg, rcolor;
//	remap(left, rleftimg, rmap[0][0], rmap[0][1], INTER_LINEAR);
//	remap(right, rrightimg, rmap[1][0], rmap[1][1], INTER_LINEAR);
//
//	////// ����У����ͼƬ
//	////imwrite("rleft.jpg", rleftimg);
//	////imwrite("rright.jpg", rrightimg);
//
//	//// �������߽��жԱ�
//	//Mat canvas;
//	//canvas.create(img_Height, img_Width * 2, CV_8UC3);
//	//hconcat(rleftimg, rrightimg, canvas);
//	//for (int j = 0; j < canvas.rows; j += 16)
//	//	line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
//	//line(canvas, Point(1280, 0), Point(1280, 720), Scalar(255, 255, 0), 1, 8);
//
//	//// ��ʾУ�����
//	//resize(canvas, canvas, Size(), 0.5, 0.5);
//	//namedWindow("rectifyResult", WINDOW_AUTOSIZE);
//	//imshow("rectifyResult", canvas);
//	//waitKey(0);
//
//
//	// �ҶȻ�
//	Mat left_img_gray, right_img_gray;
//	cvtColor(rleftimg, left_img_gray, COLOR_BGR2GRAY);
//	cvtColor(rrightimg, right_img_gray, COLOR_BGR2GRAY);
//	//// ��
//	//Mat left_img_sharpen, right_img_sharpen;
//	//Mat kernel(3, 3, CV_32F, cv::Scalar(0));
//	//kernel.at<float>(1, 1) = 5.0;
//	//kernel.at<float>(0, 1) = -1.0;
//	//kernel.at<float>(1, 0) = -1.0;
//	//kernel.at<float>(1, 2) = -1.0;
//	//kernel.at<float>(2, 1) = -1.0;
//	//filter2D(left_img_gray, left_img_sharpen, left_img_gray.depth(), kernel);
//	//filter2D(right_img_gray, right_img_sharpen, right_img_gray.depth(), kernel);
//
//	// YOLOv3 Ŀ����
//	Rect rect_left_tar(550, 400, 337, 201);
//	Rect rect_right_tar(300, 400, 337, 214);
//	Mat left_tar_img = left_img_gray(rect_left_tar); // ��ȡĿ������
//	Mat right_tar_img = right_img_gray(rect_right_tar);
//
//	// //����Ŀ������
//	//imshow("left_lines", left_tar_img);
//	//imshow("right_lines", right_tar_img);
//	////imwrite("roi1.jpg", left_tar_img);
//	////imwrite("roi2.jpg", right_tar_img);
//	//waitKey(0);
//
//	// ��Ŀ���������ֱ�߼��
//	int    length_threshold = 20; // ���ڴ���ֵ��ֱ�߱�����
//	float  distance_threshold = 1.41421356f;  // ������߶θ�Զ�ĵ㱻��Ϊ����Ⱥ��
//	double canny_th1 = 40.0; // �ݶ�С�ڴ���ֵ����Ϊ�Ǳ�Ե
//	double canny_th2 = 100.0; // �ݶȴ��ڴ���ֵһ���Ǳ�Ե
//	int    canny_aperture_size = 3;
//	bool   do_merge = true;
//	Ptr<FastLineDetector> fld = createFastLineDetector(
//		length_threshold,
//		distance_threshold,
//		canny_th1,
//		canny_th2,
//		canny_aperture_size,
//		do_merge);
//
//	vector<Vec4f> left_lines_fld;
//	left_lines_fld.push_back(Vec4f(150, 9, 186, 8));
//	left_lines_fld.push_back(Vec4f(145, 18, 153, 123));
//	left_lines_fld.push_back(Vec4f(195, 15, 212, 122));
//	left_lines_fld.push_back(Vec4f(146, 128, 91, 129));
//	left_lines_fld.push_back(Vec4f(213, 126, 277, 124));
//	left_lines_fld.push_back(Vec4f(80, 138, 81, 176));
//	left_lines_fld.push_back(Vec4f(288, 133, 294, 171));
//	left_lines_fld.push_back(Vec4f(92, 188, 285, 182));
//
//	vector<Vec4f> right_lines_fld;
//	right_lines_fld.push_back(Vec4f(86, 135, 75, 179));
//	right_lines_fld.push_back(Vec4f(82, 189, 281, 182));
//	right_lines_fld.push_back(Vec4f(288, 174, 291, 133));
//	right_lines_fld.push_back(Vec4f(189, 9, 228, 8));
//	right_lines_fld.push_back(Vec4f(182, 16, 162, 122));
//	right_lines_fld.push_back(Vec4f(233, 17, 219, 120));
//	right_lines_fld.push_back(Vec4f(98, 129, 153, 129));
//	right_lines_fld.push_back(Vec4f(223, 126, 284, 125));
//
//	// ֱ�߸�ʽת��
//	vector<Line2D> vec_left_tar_line2d;
//	vector<Line2D> vec_right_tar_line2d;
//	for (const auto& line_4f : left_lines_fld)
//	{
//		static int serial_number = 0;
//		vec_left_tar_line2d.push_back(Line2D(serial_number, line_4f));
//		++serial_number;
//	}
//	for (const auto& line_4f : right_lines_fld)
//	{
//		static int serial_number = 0;
//		vec_right_tar_line2d.push_back(Line2D(serial_number, line_4f));
//		++serial_number;
//	}
//
//	//// ����ֱ�ߺ�������
//	//Mat left_line_img = left_tar_img.clone();
//	//Mat right_line_img = right_tar_img.clone();
//	//cvtColor(left_line_img, left_line_img, COLOR_GRAY2BGR);
//	//cvtColor(right_line_img, right_line_img, COLOR_GRAY2BGR);
//	//for (const auto& line2d : vec_left_tar_line2d)
//	//{
//	//	cv::line(left_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(left_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_tar_line2d)
//	//{
//	//	cv::line(right_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(right_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//imshow("left_lines", left_line_img);
//	//imshow("right_lines", right_line_img);
//	//waitKey(0);
//
//
//	// ��ֱ�߶˵�ȡ����
//	vector<int> left_line_points_index;
//	vector<int> right_line_points_index;
//	vector<LinePoint2d> vec_left_line_points; // ������ֱ�߶˵�
//	vector<LinePoint2d> vec_right_line_points; // ������ֱ�߶˵�
//	for (const auto& line2d : vec_left_tar_line2d) // ��
//	{
//		static int index = 0;
//
//		LinePoint2d line_point1(line2d.start_point, index, line2d.serial_number, true);
//		LinePoint2d line_point2(line2d.end_point, index + 1, line2d.serial_number, false);
//
//		vec_left_line_points.push_back(line_point1);
//		vec_left_line_points.push_back(line_point2);
//		left_line_points_index.push_back(index);
//		left_line_points_index.push_back(index + 1);
//
//		index += 2;
//	}
//	for (const auto& line2d : vec_right_tar_line2d) // ��
//	{
//		static int index = 0;
//
//		LinePoint2d line_point1(line2d.start_point, index, line2d.serial_number, true);
//		LinePoint2d line_point2(line2d.end_point, index + 1, line2d.serial_number, false);
//
//		vec_right_line_points.push_back(line_point1);
//		vec_right_line_points.push_back(line_point2);
//		right_line_points_index.push_back(index);
//		right_line_points_index.push_back(index + 1);
//
//		index += 2;
//	}
//
//	// �Զ˵������ֱ�߽��з���
//	vector<vector<int>> vec_left_near_points_group; // �洢����vec_left_line_points���±�
//	vector<vector<int>> vec_right_near_points_group;
//
//	double line_point_near_threshold = 25; // �˵��������ֵ
//
//	while (left_line_points_index.size() != 0)
//	{
//		int curr_index = left_line_points_index[0];
//		vector<int> curr_near;
//		LinePoint2d curr_point = vec_left_line_points[curr_index];
//		for (int j = 0; j < left_line_points_index.size(); j++)
//		{
//			int comp_index = left_line_points_index[j];
//
//			LinePoint2d comp_point = vec_left_line_points[comp_index];
//
//			double dis = computeTowPointDistance(curr_point.point_2d, comp_point.point_2d);
//
//			if (dis <= line_point_near_threshold)
//			{
//				curr_near.push_back(comp_index);
//			}
//		}
//
//		// ɾ���Ѿ���ԵĶ˵�
//		for (int i = 0; i < curr_near.size(); i++)
//		{
//			deleteElement(left_line_points_index, curr_near[i]);
//		}
//		vec_left_near_points_group.push_back(curr_near);
//	}
//
//	while (right_line_points_index.size() != 0)
//	{
//		int curr_index = right_line_points_index[0];
//		vector<int> curr_near;
//		LinePoint2d curr_point = vec_right_line_points[curr_index];
//		for (int j = 0; j < right_line_points_index.size(); j++)
//		{
//			int comp_index = right_line_points_index[j];
//
//			LinePoint2d comp_point = vec_right_line_points[comp_index];
//
//			double dis = computeTowPointDistance(curr_point.point_2d, comp_point.point_2d);
//
//			if (dis <= line_point_near_threshold)
//			{
//				curr_near.push_back(comp_index);
//			}
//		}
//
//		// ɾ���Ѿ���ԵĶ˵�
//		for (int i = 0; i < curr_near.size(); i++)
//		{
//			deleteElement(right_line_points_index, curr_near[i]);
//		}
//		vec_right_near_points_group.push_back(curr_near);
//	}
//
//	// ֱ�ߺϲ�
//	double line_slope_similarity_threshold = 0.5; // ֱ��б�����ƶ���ֵ
//
//	// �Զ˵������ֱ�߽����ӳ����󽻵�
//	Mat left_cross_point_image = left_tar_img.clone();
//	Mat right_cross_point_image = right_tar_img.clone();
//	cvtColor(left_cross_point_image, left_cross_point_image, COLOR_GRAY2BGR);
//	cvtColor(right_cross_point_image, right_cross_point_image, COLOR_GRAY2BGR);
//	for (const auto& near_points : vec_left_near_points_group) // ��
//	{
//		if (near_points.size() == 1)
//			continue;
//		if (near_points.size() == 2)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			LinePoint2d point1 = vec_left_line_points[point_index1];
//			LinePoint2d point2 = vec_left_line_points[point_index2];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			Line2D& line1 = vec_left_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_left_tar_line2d[line_serial_number2];
//			Point2f cross_point = conputeIntersectionPoint(line1, line2);
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			cv::circle(left_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//		if (near_points.size() == 3)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			int point_index3 = near_points[2];
//			LinePoint2d point1 = vec_left_line_points[point_index1];
//			LinePoint2d point2 = vec_left_line_points[point_index2];
//			LinePoint2d point3 = vec_left_line_points[point_index3];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			int line_serial_number3 = point3.line_serial_number;
//			Line2D& line1 = vec_left_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_left_tar_line2d[line_serial_number2];
//			Line2D& line3 = vec_left_tar_line2d[line_serial_number3];
//			Point2f cross_point1 = conputeIntersectionPoint(line1, line2);
//			Point2f cross_point2 = conputeIntersectionPoint(line1, line3);
//			Point2f cross_point3 = conputeIntersectionPoint(line2, line3);
//			Point2f cross_point = (cross_point1 + cross_point2 + cross_point3) / 3.0;
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			if (point3.is_start_point == true)
//				line3.start_point = cross_point;
//			else
//				line3.end_point = cross_point;
//
//			cv::circle(left_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//	}
//	for (const auto& near_points : vec_right_near_points_group) // ��
//	{
//		if (near_points.size() == 1)
//			continue;
//		if (near_points.size() == 2)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			LinePoint2d point1 = vec_right_line_points[point_index1];
//			LinePoint2d point2 = vec_right_line_points[point_index2];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			Line2D& line1 = vec_right_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_right_tar_line2d[line_serial_number2];
//			Point2f cross_point = conputeIntersectionPoint(line1, line2);
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			cv::circle(right_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//		if (near_points.size() == 3)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			int point_index3 = near_points[2];
//			LinePoint2d point1 = vec_right_line_points[point_index1];
//			LinePoint2d point2 = vec_right_line_points[point_index2];
//			LinePoint2d point3 = vec_right_line_points[point_index3];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			int line_serial_number3 = point3.line_serial_number;
//			Line2D& line1 = vec_right_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_right_tar_line2d[line_serial_number2];
//			Line2D& line3 = vec_right_tar_line2d[line_serial_number3];
//			Point2f cross_point1 = conputeIntersectionPoint(line1, line2);
//			Point2f cross_point2 = conputeIntersectionPoint(line1, line3);
//			Point2f cross_point3 = conputeIntersectionPoint(line2, line3);
//			Point2f cross_point = (cross_point1 + cross_point2 + cross_point3) / 3.0;
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			if (point3.is_start_point == true)
//				line3.start_point = cross_point;
//			else
//				line3.end_point = cross_point;
//
//			cv::circle(right_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//	}
//
//	//// ��ʾֱ�߽���
//	//imshow("left_cross_point_image", left_cross_point_image);
//	//imshow("right_cross_point_image", right_cross_point_image);
//	//waitKey(0);
//
//	//// ��ʾ����������Ϊ�˵��ֱ��
//	//Mat left_line_img = left_tar_img.clone();
//	//Mat right_line_img = right_tar_img.clone();
//	//cvtColor(left_line_img, left_line_img, COLOR_GRAY2BGR);
//	//cvtColor(right_line_img, right_line_img, COLOR_GRAY2BGR);
//	//for (const auto& line2d : vec_left_tar_line2d)
//	//{
//	//	cv::line(left_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(left_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_tar_line2d)
//	//{
//	//	cv::line(right_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(right_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//imshow("left_lines", left_line_img);
//	//imshow("right_lines", right_line_img);
//	//waitKey(0);
//
//
//	// ֱ�߾ֲ�����ӳ�䵽ȫ�����꣬���¼���б��
//	vector<Line2D> vec_left_line2d;
//	vector<Line2D> vec_right_line2d;
//	for (const auto& line2d : vec_left_tar_line2d) // ��
//	{
//		static int index = 0;
//		float x1 = line2d.start_point.x + rect_left_tar.x;
//		float y1 = line2d.start_point.y + rect_left_tar.y;
//		float x2 = line2d.end_point.x + rect_left_tar.x;
//		float y2 = line2d.end_point.y + rect_left_tar.y;
//		Vec4f line4f(x1, y1, x2, y2);
//		vec_left_line2d.push_back(Line2D(index, line4f));
//		++index;
//	}
//	for (const auto& line2d : vec_right_tar_line2d) // ��
//	{
//		static int index = 0;
//		float x1 = line2d.start_point.x + rect_right_tar.x;
//		float y1 = line2d.start_point.y + rect_right_tar.y;
//		float x2 = line2d.end_point.x + rect_right_tar.x;
//		float y2 = line2d.end_point.y + rect_right_tar.y;
//		Vec4f line4f(x1, y1, x2, y2);
//		vec_right_line2d.push_back(Line2D(index, line4f));
//		++index;
//	}
//
//	//// ���ƴ���õ�ȫ��ֱ��
//	//Mat left_global_line_img = rleftimg.clone();
//	//Mat right_global_line_img = rrightimg.clone();
//
//	//for (const auto& line2d : vec_left_line2d)
//	//{
//	//	cv::line(left_global_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 2);
//	//	cv::putText(left_global_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_line2d)
//	//{
//	//	cv::line(right_global_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 2);
//	//	cv::putText(right_global_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//
//	//Mat left_right_global_line_img; // ƴ������ͼ��
//	//left_right_global_line_img.create(img_Height, img_Width * 2, CV_8UC3);
//	//hconcat(left_global_line_img, right_global_line_img, left_right_global_line_img);
//	//imshow("left_global_line_img", left_global_line_img);
//	//imshow("right_global_line_img", right_global_line_img);
//	//imshow("left_right_global_line_img", left_right_global_line_img);
//	//waitKey(0);
//
//	// ֱ��ƥ��
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_reconstruct(new pcl::PointCloud<pcl::PointXYZ>()); // ֱ�ߵ���
//	std::vector<std::vector<int>> vec_match_index_group; // ƥ������
//	//lineMatch(vec_left_line2d, vec_right_line2d, vec_match_index_group);
//	std::vector<int> match1; match1.push_back(0); match1.push_back(3); vec_match_index_group.push_back(match1);
//	std::vector<int> match2; match2.push_back(1); match2.push_back(4); vec_match_index_group.push_back(match2);
//	std::vector<int> match3; match3.push_back(2); match3.push_back(5); vec_match_index_group.push_back(match3);
//	//std::vector<int> match4; match4.push_back(3); match4.push_back(6); vec_match_index_group.push_back(match4);/////////
//	std::vector<int> match5; match5.push_back(4); match5.push_back(7); vec_match_index_group.push_back(match5);
//	std::vector<int> match6; match6.push_back(5); match6.push_back(0); vec_match_index_group.push_back(match6);
//	std::vector<int> match7; match7.push_back(6); match7.push_back(2); vec_match_index_group.push_back(match7);
//	std::vector<int> match8; match8.push_back(7); match8.push_back(1); vec_match_index_group.push_back(match8);
//
//	// ֱ���ؽ�
//	vector<Line3D> vec_line3d;
//	for (const auto vec_match_index : vec_match_index_group)
//	{
//		if (vec_match_index.size() != 2)
//			continue;
//
//		int left_index = vec_match_index[0];
//		int right_index = vec_match_index[1];
//		const Line2D& left_line = vec_left_line2d[left_index];
//		const Line2D& right_line = vec_right_line2d[right_index];
//
//		// p1��㣬p2�յ�
//		const Point2f& left_p1 = left_line.start_point;
//		const Point2f& left_p2 = left_line.end_point;
//		const Point2f& right_p1 = right_line.start_point;
//		const Point2f& right_p2 = right_line.end_point;
//
//		// ����ռ�ֱ�������˵�
//		const float f = 622.812; // ���࣬���㣬���߾���
//		const float cx = 654.193;
//		const float cy = 344.569;
//		const float Tx = 60.14;
//
//		float p1_disp = left_p1.x - right_p1.x;
//
//		float p1_z_w = f * Tx / p1_disp;
//		float p1_x_w = (left_p1.x - cx)*p1_z_w / f;
//		float p1_y_w = (left_p1.y - cy)*p1_z_w / f;
//
//		float p2_disp = left_p2.x - right_p2.x;
//		float p2_z_w = f * Tx / p2_disp;
//		float p2_x_w = (left_p2.x - cx)*p2_z_w / f;
//		float p2_y_w = (left_p2.y - cy)*p2_z_w / f;
//
//		// push�ؽ��Ŀռ�ֱ��
//		Line3D line_3d(Point3f(p1_x_w, p1_y_w, p1_z_w), Point3f(p2_x_w, p2_y_w, p2_z_w));
//		vec_line3d.push_back(line_3d);
//
//		// push �ռ�ֱ�߶˵�
//		pcl::PointXYZ p1_w(p1_x_w, p1_y_w, p1_z_w);
//		pcl::PointXYZ p2_w(p2_x_w, p2_y_w, p2_z_w);
//		cloud_reconstruct->points.push_back(p1_w);
//		cloud_reconstruct->points.push_back(p2_w);
//	}
//
//	// push δ�ؽ��ĵ�
//	Line3D line_3d(Point3f(12.186, 45.412, 154.33), Point3f(-6.06, 45.833, 154.6));
//	vec_line3d.push_back(line_3d);
//
//	// ȫ��ֱ�߲����ɵ���
//	float sample_dis = 1.0; // ֱ�ߴ�Ͼ���
//
//	for (const auto& line_3d : vec_line3d)
//	{
//		// ֱ�ߵȾ�����
//		float x0 = line_3d.x0_vec.x;
//		float y0 = line_3d.x0_vec.y;
//		float z0 = line_3d.x0_vec.z;
//
//		float v1 = line_3d.normal_vec.x;
//		float v2 = line_3d.normal_vec.y;
//		float v3 = line_3d.normal_vec.z;
//
//		float v_square = line_3d.normal_vec.dot(line_3d.normal_vec);
//		float deta_t = sample_dis / std::sqrt(v_square);
//
//		for (float t = line_3d.t_start; t < line_3d.t_end; t += deta_t)
//		{
//			pcl::PointXYZ point_3d;
//			point_3d.x = x0 + v1 * t;
//			point_3d.y = y0 + v2 * t;
//			point_3d.z = z0 + v3 * t;
//			cloud_reconstruct->points.push_back(point_3d);
//		}
//	}
//
//	// ������׼
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src(new pcl::PointCloud<pcl::PointXYZ>());
//	createTWireFramePointCloud(cloud_src);
//	// ������Ʊ��淨��
//	clock_t fpfh_start_time = clock();
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_src; // Դ
//	ne_src.setInputCloud(cloud_src);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_src.setSearchMethod(tree_src);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_src.setRadiusSearch(2.0);
//	ne_src.setNumberOfThreads(8);
//	ne_src.compute(*cloud_src_normals);
//
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_tar; // �ؽ�
//	ne_tar.setInputCloud(cloud_reconstruct);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_reconstruct(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_tar.setSearchMethod(tree_reconstruct);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_reconstruct_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_tar.setRadiusSearch(2.0);
//	ne_src.setNumberOfThreads(8);
//	ne_tar.compute(*cloud_reconstruct_normals);
//
//	//����FPFH
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_src; // Դ
//	fpfh_src.setInputCloud(cloud_src);
//	fpfh_src.setInputNormals(cloud_src_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_src_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_src.setSearchMethod(tree_src_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_src(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_src.setRadiusSearch(1.8);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_src.compute(*fpfhs_src);
//	std::cout << "compute *cloud_src fpfh" << endl;
//
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_tar; // �ؽ�
//	fpfh_tar.setInputCloud(cloud_reconstruct);
//	fpfh_tar.setInputNormals(cloud_reconstruct_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_reconstruct_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_tar.setSearchMethod(tree_reconstruct_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_reconstruct(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_tar.setRadiusSearch(1.8);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_tar.compute(*fpfhs_reconstruct);
//
//	std::cout << "compute *cloud_tar fpfh" << endl;
//	clock_t fpfh_end_time = clock();
//	cout << "fpfh time: " << (double)(fpfh_end_time - fpfh_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//SAC��׼
//	clock_t sac_start_time = clock();
//
//	pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> scia;
//	scia.setInputSource(cloud_src);
//	scia.setInputTarget(cloud_reconstruct);
//	scia.setSourceFeatures(fpfhs_src);
//	scia.setTargetFeatures(fpfhs_reconstruct);
//
//	//scia.setMinSampleDistance(1);
//	//scia.setNumberOfSamples(2);
//	//scia.setCorrespondenceRandomness(20);
//
//	pcl::PointCloud<pcl::PointXYZ>::Ptr sac_result(new pcl::PointCloud<pcl::PointXYZ>);
//	scia.align(*sac_result);
//	Eigen::Matrix4f sac_trans = scia.getFinalTransformation();
//	clock_t sac_end_time = clock();
//	cout << "sac time: " << (double)(sac_end_time - sac_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//icp��׼
//	clock_t icp_start_time = clock();
//	pcl::PointCloud<pcl::PointXYZ>::Ptr icp_result(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
//	icp.setInputSource(cloud_src);
//	icp.setInputTarget(cloud_reconstruct);
//	// ��������Ӧ����� (e.g., correspondences with higher distances will be ignored)
//	icp.setMaxCorrespondenceDistance(40);
//	// ����������
//	icp.setMaximumIterations(100);
//	// ���α仯����֮��Ĳ�ֵ
//	icp.setTransformationEpsilon(1e-10);
//	// �������
//	icp.setEuclideanFitnessEpsilon(0.2);
//	icp.align(*icp_result, sac_trans);
//	clock_t icp_end_time = clock();
//	cout << "icp time: " << (double)(icp_end_time - icp_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//	// ��ȡ���յı任����
//	Eigen::Matrix4f final_trans_pcl = icp.getFinalTransformation();
//	cv::Mat_<float> final_trans_cv(4, 4);
//	for (size_t i = 0; i < 4; i++)
//	{
//		for (size_t j = 0; j < 4; j++)
//		{
//			final_trans_cv(i, j) = final_trans_pcl(i, j);
//		}
//	}
//
//	// ��������ϵ������
//	vector<cv::Point3f> vec_axis_3d;
//	vector<cv::Point2f> vec_axis_2d;
//	vec_axis_3d.push_back(Point3f(10.0, 0.0, 0.0));
//	vec_axis_3d.push_back(Point3f(0.0, 10.0, 0.0));
//	vec_axis_3d.push_back(Point3f(0.0, 0.0, 10.0));
//	vec_axis_3d.push_back(Point3f(0.0, 0.0, 0.0));
//
//	// ��������ϵת��������
//	Eigen::Vector3f final_angle;
//	Eigen::Vector3f final_translate;
//	matrix2Angle(final_trans_pcl, final_angle);
//	matrix2Translate(final_trans_pcl, final_translate);
//	cv::Mat rVec(3, 1, cv::DataType<float>::type); // Rotation vector
//	rVec.at<float>(0) = final_angle.x();
//	rVec.at<float>(1) = final_angle.y();
//	rVec.at<float>(2) = final_angle.z();
//
//	cv::Mat tVec(3, 1, cv::DataType<float>::type); // Translation vector
//	tVec.at<float>(0) = final_translate.x();
//	tVec.at<float>(1) = final_translate.y();
//	tVec.at<float>(2) = final_translate.z();
//
//	cout << "λ�ƣ�" << tVec << endl;
//	cout << "��̬��" << rVec << endl;
//
//	cv::projectPoints(vec_axis_3d, rVec, tVec, Pl.colRange(0, 3), Mat(), vec_axis_2d);
//
//	// ������ͶӰ���������
//	Mat reproject_axis_img = rleftimg.clone();
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[0], Scalar(0, 0, 255), 2);
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[1], Scalar(0, 255, 0), 2);
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[2], Scalar(255, 0, 0), 2);
//	imshow("reproject_axis_img", reproject_axis_img);
//	waitKey(0);
//
//	// ��ʾ����
//	// ��ɫ��Դ����
//	// ��ɫ��Ŀ�����
//	// ��ɫ����׼����
//	pcl::visualization::PCLVisualizer viewer("registration Viewer");
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h(cloud_src, 0, 255, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> tgt_h(cloud_reconstruct, 255, 0, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> final_h(icp_result, 0, 0, 255);
//	viewer.setBackgroundColor(1.0, 1.0, 1.0);
//
//	//viewer.addPointCloud(cloud_src, src_h, "tar_points");
//	viewer.addPointCloud(cloud_reconstruct, tgt_h, "src_points");
//	viewer.addPointCloud(icp_result, final_h, "icp_points");
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "src_points"); // ���õ��ƴ�С
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "tar_points");
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "icp_points");
//	viewer.addCoordinateSystem(20.0); // �������ָʾ�ᣬ���ű���-------��ɫX����ɫY����ɫZ
//
//	// �����ؽ�����
//	cloud_reconstruct->width = 1;
//	cloud_reconstruct->height = cloud_reconstruct->points.size();
//	pcl::io::savePCDFileASCII("./cloud_reconstruct_src/t_reconstruct.pcd", *cloud_reconstruct);
//	cout << "cloud_reconstruct is saved" << endl;
//
//	// ����Դ����
//	cloud_src->width = 1;
//	cloud_src->height = cloud_src->points.size();
//	pcl::io::savePCDFileASCII("./cloud_reconstruct_src/t_src.pcd", *cloud_src);
//	cout << "cloud_src is saved" << endl;
//
//	// ������׼����
//	icp_result->width = 1;
//	icp_result->height = icp_result->points.size();
//	pcl::io::savePCDFileASCII("./cloud_reconstruct_src/t_icp_result.pcd", *icp_result);
//	cout << "icp_result is saved" << endl;
//
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(100);
//		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	}
//
//	return 0;
//}
//
//
//// T����
//void createTWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame)
//{
//
//	// 50
//	for (int i = 0; i <= 50; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = -25 + i; p1.y = 0; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//	}
//
//
//	// 18
//	for (int i = 0; i <= 18; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = -25 + i; p1.y = 14; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 7 + i; p2.y = 14; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//	}
//
//	// 14
//	for (int i = 0; i <= 14; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = -25; p1.y = i; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 25; p2.y = i; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = i - 7; p3.y = 50; p3.z = 0;
//		cloud_wire_frame->points.push_back(p3);
//
//	}
//
//	// 36
//	for (int i = 0; i <= 36; i++)
//	{
//		pcl::PointXYZ p1;
//		p1.x = -7; p1.y = 14 + i; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 7; p2.y = 14 + i; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//	}
//}
//
//
////����תƽ�ƾ��������ת�Ƕ�
//void matrix2Angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle)
//{
//	double ax, ay, az;
//	if (result_trans(2, 0) == 1 || result_trans(2, 0) == -1)
//	{
//		az = 0;
//		double dlta;
//		dlta = atan2(result_trans(0, 1), result_trans(0, 2));
//		if (result_trans(2, 0) == -1)
//		{
//			ay = M_PI / 2;
//			ax = az + dlta;
//		}
//		else
//		{
//			ay = -M_PI / 2;
//			ax = -az + dlta;
//		}
//	}
//	else
//	{
//		ay = -asin(result_trans(2, 0));
//		ax = atan2(result_trans(2, 1) / cos(ay), result_trans(2, 2) / cos(ay));
//		az = atan2(result_trans(1, 0) / cos(ay), result_trans(0, 0) / cos(ay));
//	}
//	result_angle << ax, ay, az;
//}
//
////����תƽ�ƾ������ƽ�ƾ���
//void matrix2Translate(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_translate)
//{
//	result_translate << result_trans(0, 3), result_trans(1, 3), result_trans(2, 3);
//}


///*
//  ������
//*/
//
//#pragma once
//#include "LineAlgorithm.h"
//
//#include <iostream>  
//#include <time.h>
//
//#include <opencv2/opencv.hpp>
//#include <opencv2/ximgproc.hpp>
//#include <opencv2/xfeatures2d.hpp>
//
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/common/io.h>
//#include <pcl/point_cloud.h>
//
//#include <pcl/PolygonMesh.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <pcl/io/vtk_lib_io.h>
//#include <pcl/io/io.h>
//#include <pcl/io/pcd_io.h>//pcd ��д����ص�ͷ�ļ���
//#include <pcl/point_types.h> //PCL��֧�ֵĵ�����ͷ�ļ�
//#include <pcl/registration/ia_ransac.h>
//#include <pcl/point_cloud.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/features/normal_3d_omp.h>
//#include <pcl/features/fpfh.h>
//#include <pcl/features/fpfh_omp.h>
//#include <pcl/search/kdtree.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/filter.h>
//#include <pcl/registration/icp.h>
//#include <pcl/visualization/pcl_visualizer.h>
//
//using namespace cv;
//using namespace cv::ximgproc;
//using namespace std;
//
//// ͼƬ��С
//const int img_Width = 1280;
//const int img_Height = 720;
//
////У����ת����R��ͶӰ����P ��ͶӰ����Q
//Mat Rl, Rr, Pl, Pr, Q;
//
//// ����ģ�����
//void createCuboidWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame);
//
//// ����ת��תƽ������
//void matrix2Angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle);
//void matrix2Translate(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_translate);
//
//int main()
//{
//
//	Mat left_and_right = imread("./left_right_img/IMG_0022.jpg");
//	Rect rect_left(0, 0, img_Width, img_Height);
//	Rect rect_right(img_Width, 0, img_Width, img_Height);
//
//	Mat left = left_and_right(rect_left);
//	Mat right = left_and_right(rect_right);
//
//	//imshow("src_left", left);
//	//imshow("src_right", right);
//	//waitKey(0);
//
//	Mat cameraMatrix[2], distCoeffs[2];
//	Mat Q;
//
//	// ����ڲ�
//	cameraMatrix[0] = Mat((Mat_<double>(3, 3) << 731.5022, 0.0, 0.0,
//		0.0, 732.0644, 0.0,
//		674.6773, 357.5894, 1.0)).t();
//
//	cameraMatrix[1] = Mat((Mat_<double>(3, 3) << 727.3993, 0.0, 0.0,
//		0.0, 727.6805, 0.0,
//		634.0939, 339.6110, 1.0)).t();
//
//	// ����ϵ��
//	distCoeffs[0] = (Mat_<double>(4, 1) << 0.1070, -0.1208, 0.0, -0.0004);
//	distCoeffs[1] = (Mat_<double>(4, 1) << 0.1145, -0.1517, 0.0007, -0.0004);
//
//	// ������λ��
//	Mat R = Mat((Mat_<double>(3, 3) << 1.0, 4.7372e-04, 4.8307e-04,
//		-4.7146e-04, 1.0, -0.0047,
//		-0.0005, 0.0047, 1.0)).t();
//	Mat T = (Mat_<double>(3, 1) << -60.1450, -0.0701, -0.0956);
//
//	// ͼ��У��֮�����Ч����
//	Rect validRoi[2];
//
//	// bougust����У��
//	stereoRectify(cameraMatrix[0], distCoeffs[0],
//		cameraMatrix[1], distCoeffs[1],
//		Size(img_Width, img_Height), R, T, Rl, Rr, Pl, Pr, Q,
//		CALIB_ZERO_DISPARITY, 1, Size(img_Width, img_Height), &validRoi[0], &validRoi[1]); //validPixROI1  validPixROI2 - ��ѡ�����������Rect�����ݡ����ڲ����������ض���Ч
//
//	//��ʾ���������ͶӰ����
//	cout << "camL" << cameraMatrix[0] << endl;
//	cout << "camR" << cameraMatrix[1] << endl;
//	cout << "Pro1" << Pl << endl;
//	cout << "Pro2" << Pr << endl;
//	cout << "Q" << Q << endl;
//
//	// ����������ͼ��У������ӳ���
//	Mat rmap[2][2];
//	initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], Rl, Pl, Size(img_Width, img_Height), CV_16SC2, rmap[0][0], rmap[0][1]);
//	initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], Rr, Pr, Size(img_Width, img_Height), CV_16SC2, rmap[1][0], rmap[1][1]);
//
//	// ��Դͼ�������ӳ�䣬ӳ�����ΪУ�����ͼ��
//	Mat rleftimg, rrightimg, rcolor;
//	remap(left, rleftimg, rmap[0][0], rmap[0][1], INTER_LINEAR);
//	remap(right, rrightimg, rmap[1][0], rmap[1][1], INTER_LINEAR);
//
//	////// ����У����ͼƬ
//	////imwrite("rleft.jpg", rleftimg);
//	////imwrite("rright.jpg", rrightimg);
//
//	//// �������߽��жԱ�
//	//Mat canvas;
//	//canvas.create(img_Height, img_Width * 2, CV_8UC3);
//	//hconcat(rleftimg, rrightimg, canvas);
//	//for (int j = 0; j < canvas.rows; j += 16)
//	//	line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
//	//line(canvas, Point(1280, 0), Point(1280, 720), Scalar(255, 255, 0), 1, 8);
//
//	//// ��ʾУ�����
//	//resize(canvas, canvas, Size(), 0.5, 0.5);
//	//namedWindow("rectifyResult", WINDOW_AUTOSIZE);
//	//imshow("rectifyResult", canvas);
//	//waitKey(0);
//
//
//	// �ҶȻ�
//	Mat left_img_gray, right_img_gray;
//	cvtColor(rleftimg, left_img_gray, COLOR_BGR2GRAY);
//	cvtColor(rrightimg, right_img_gray, COLOR_BGR2GRAY);
//	//// ��
//	//Mat left_img_sharpen, right_img_sharpen;
//	//Mat kernel(3, 3, CV_32F, cv::Scalar(0));
//	//kernel.at<float>(1, 1) = 5.0;
//	//kernel.at<float>(0, 1) = -1.0;
//	//kernel.at<float>(1, 0) = -1.0;
//	//kernel.at<float>(1, 2) = -1.0;
//	//kernel.at<float>(2, 1) = -1.0;
//	//filter2D(left_img_gray, left_img_sharpen, left_img_gray.depth(), kernel);
//	//filter2D(right_img_gray, right_img_sharpen, right_img_gray.depth(), kernel);
//
//	// YOLOv3 Ŀ����
//	Rect rect_left_tar(305, 180, 337, 230);
//	Rect rect_right_tar(120, 180, 337, 230);
//	Mat left_tar_img = left_img_gray(rect_left_tar); // ��ȡĿ������
//	Mat right_tar_img = right_img_gray(rect_right_tar);
//
//	////����Ŀ������
//	//imshow("left_lines", left_tar_img);
//	//imshow("right_lines", right_tar_img);
//	////imwrite("roi1.jpg", left_tar_img);
//	////imwrite("roi2.jpg", right_tar_img);
//	//waitKey(0);
//
//	// ��Ŀ���������ֱ�߼��
//	int    length_threshold = 20; // ���ڴ���ֵ��ֱ�߱�����
//	float  distance_threshold = 1.41421356f;  // ������߶θ�Զ�ĵ㱻��Ϊ����Ⱥ��
//	double canny_th1 = 40.0; // �ݶ�С�ڴ���ֵ����Ϊ�Ǳ�Ե
//	double canny_th2 = 100.0; // �ݶȴ��ڴ���ֵһ���Ǳ�Ե
//	int    canny_aperture_size = 3;
//	bool   do_merge = true;
//	Ptr<FastLineDetector> fld = createFastLineDetector(
//		length_threshold,
//		distance_threshold,
//		canny_th1,
//		canny_th2,
//		canny_aperture_size,
//		do_merge);
//
//
//	//vector<Vec4f> left_lines_fld;
//	//left_lines_fld.push_back(Vec4f(150, 9, 186, 8));
//	//left_lines_fld.push_back(Vec4f(145, 18, 153, 123));
//	//left_lines_fld.push_back(Vec4f(195, 15, 212, 122));
//	//left_lines_fld.push_back(Vec4f(146, 128, 91, 129));
//	//left_lines_fld.push_back(Vec4f(213, 126, 277, 124));
//	//left_lines_fld.push_back(Vec4f(80, 138, 81, 176));
//	//left_lines_fld.push_back(Vec4f(288, 133, 294, 171));
//	//left_lines_fld.push_back(Vec4f(92, 188, 285, 182));
//
//	//vector<Vec4f> right_lines_fld;
//	//right_lines_fld.push_back(Vec4f(86, 135, 75, 179));
//	//right_lines_fld.push_back(Vec4f(82, 189, 281, 182));
//	//right_lines_fld.push_back(Vec4f(288, 174, 291, 133));
//	//right_lines_fld.push_back(Vec4f(189, 9, 228, 8));
//	//right_lines_fld.push_back(Vec4f(182, 16, 162, 122));
//	//right_lines_fld.push_back(Vec4f(233, 17, 219, 120));
//	//right_lines_fld.push_back(Vec4f(98, 129, 153, 129));
//	//right_lines_fld.push_back(Vec4f(223, 126, 284, 125));
//
//	vector<Vec4f> left_lines_fld;
//	vector<Vec4f> right_lines_fld;
//	fld->detect(left_tar_img, left_lines_fld);
//	fld->detect(right_tar_img, right_lines_fld);
//	// ֱ�߸�ʽת��
//	vector<Line2D> vec_left_tar_line2d;
//	vector<Line2D> vec_right_tar_line2d;
//	for (const auto& line_4f : left_lines_fld) // ��
//	{
//		static int serial_number = 0;
//		if (serial_number == 0)
//		{
//			vec_left_tar_line2d.push_back(Line2D(serial_number, Vec4f(62, 88, 175, 151)));
//			++serial_number;
//			continue;
//		}
//		if (serial_number == 6)
//		{
//			vec_left_tar_line2d.push_back(Line2D(serial_number, Vec4f(177.59, 152.719, 273, 73)));
//			++serial_number;
//			continue;
//		}
//
//		vec_left_tar_line2d.push_back(Line2D(serial_number, line_4f));
//		++serial_number;
//	}
//
//	for (const auto& line_4f : right_lines_fld) // ��
//	{
//		static int serial_number = 0;
//
//		if (serial_number == 0)
//		{
//			vec_right_tar_line2d.push_back(Line2D(serial_number, Vec4f(254, 73, 137, 148)));
//			++serial_number;
//			continue;
//		}
//		if (serial_number == 7)
//		{
//			vec_right_tar_line2d.push_back(Line2D(serial_number, Vec4f(39, 86, 125, 148)));
//			++serial_number;
//			continue;
//		}
//
//		vec_right_tar_line2d.push_back(Line2D(serial_number, line_4f));
//		++serial_number;
//	}
//
//	vec_right_tar_line2d.push_back(Line2D(8, Vec4f(133, 151, 182, 219)));
//
//	//// ����ֱ�ߺ�������
//	//Mat left_line_img = left_tar_img.clone();
//	//Mat right_line_img = right_tar_img.clone();
//	//cvtColor(left_line_img, left_line_img, COLOR_GRAY2BGR);
//	//cvtColor(right_line_img, right_line_img, COLOR_GRAY2BGR);
//	//for (const auto& line2d : vec_left_tar_line2d)
//	//{
//	//	cv::line(left_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(left_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_tar_line2d)
//	//{
//	//	cv::line(right_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(right_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//imshow("left_lines", left_line_img);
//	//imshow("right_lines", right_line_img);
//	//waitKey(0);
//
//
//	// ��ֱ�߶˵�ȡ����
//	vector<int> left_line_points_index;
//	vector<int> right_line_points_index;
//	vector<LinePoint2d> vec_left_line_points; // ������ֱ�߶˵�
//	vector<LinePoint2d> vec_right_line_points; // ������ֱ�߶˵�
//	for (const auto& line2d : vec_left_tar_line2d) // ��
//	{
//		static int index = 0;
//
//		LinePoint2d line_point1(line2d.start_point, index, line2d.serial_number, true);
//		LinePoint2d line_point2(line2d.end_point, index + 1, line2d.serial_number, false);
//
//		vec_left_line_points.push_back(line_point1);
//		vec_left_line_points.push_back(line_point2);
//		left_line_points_index.push_back(index);
//		left_line_points_index.push_back(index + 1);
//
//		index += 2;
//	}
//	for (const auto& line2d : vec_right_tar_line2d) // ��
//	{
//		static int index = 0;
//
//		LinePoint2d line_point1(line2d.start_point, index, line2d.serial_number, true);
//		LinePoint2d line_point2(line2d.end_point, index + 1, line2d.serial_number, false);
//
//		vec_right_line_points.push_back(line_point1);
//		vec_right_line_points.push_back(line_point2);
//		right_line_points_index.push_back(index);
//		right_line_points_index.push_back(index + 1);
//
//		index += 2;
//	}
//
//	// �Զ˵������ֱ�߽��з���
//	vector<vector<int>> vec_left_near_points_group; // �洢����vec_left_line_points���±�
//	vector<vector<int>> vec_right_near_points_group;
//
//	double line_point_near_threshold = 25; // �˵��������ֵ
//
//	while (left_line_points_index.size() != 0)
//	{
//		int curr_index = left_line_points_index[0];
//		vector<int> curr_near;
//		LinePoint2d curr_point = vec_left_line_points[curr_index];
//		for (int j = 0; j < left_line_points_index.size(); j++)
//		{
//			int comp_index = left_line_points_index[j];
//
//			LinePoint2d comp_point = vec_left_line_points[comp_index];
//
//			double dis = computeTowPointDistance(curr_point.point_2d, comp_point.point_2d);
//
//			if (dis <= line_point_near_threshold)
//			{
//				curr_near.push_back(comp_index);
//			}
//		}
//
//		// ɾ���Ѿ���ԵĶ˵�
//		for (int i = 0; i < curr_near.size(); i++)
//		{
//			deleteElement(left_line_points_index, curr_near[i]);
//		}
//		vec_left_near_points_group.push_back(curr_near);
//	}
//
//	while (right_line_points_index.size() != 0)
//	{
//		int curr_index = right_line_points_index[0];
//		vector<int> curr_near;
//		LinePoint2d curr_point = vec_right_line_points[curr_index];
//		for (int j = 0; j < right_line_points_index.size(); j++)
//		{
//			int comp_index = right_line_points_index[j];
//
//			LinePoint2d comp_point = vec_right_line_points[comp_index];
//
//			double dis = computeTowPointDistance(curr_point.point_2d, comp_point.point_2d);
//
//			if (dis <= line_point_near_threshold)
//			{
//				curr_near.push_back(comp_index);
//			}
//		}
//
//		// ɾ���Ѿ���ԵĶ˵�
//		for (int i = 0; i < curr_near.size(); i++)
//		{
//			deleteElement(right_line_points_index, curr_near[i]);
//		}
//		vec_right_near_points_group.push_back(curr_near);
//	}
//
//	// ֱ�ߺϲ�
//	double line_slope_similarity_threshold = 0.5; // ֱ��б�����ƶ���ֵ
//
//	// �Զ˵������ֱ�߽����ӳ����󽻵�
//	Mat left_cross_point_image = left_tar_img.clone();
//	Mat right_cross_point_image = right_tar_img.clone();
//	cvtColor(left_cross_point_image, left_cross_point_image, COLOR_GRAY2BGR);
//	cvtColor(right_cross_point_image, right_cross_point_image, COLOR_GRAY2BGR);
//	for (const auto& near_points : vec_left_near_points_group) // ��
//	{
//		if (near_points.size() == 1)
//			continue;
//		if (near_points.size() == 2)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			LinePoint2d point1 = vec_left_line_points[point_index1];
//			LinePoint2d point2 = vec_left_line_points[point_index2];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			Line2D& line1 = vec_left_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_left_tar_line2d[line_serial_number2];
//			Point2f cross_point = conputeIntersectionPoint(line1, line2);
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			cv::circle(left_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//		if (near_points.size() == 3)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			int point_index3 = near_points[2];
//			LinePoint2d point1 = vec_left_line_points[point_index1];
//			LinePoint2d point2 = vec_left_line_points[point_index2];
//			LinePoint2d point3 = vec_left_line_points[point_index3];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			int line_serial_number3 = point3.line_serial_number;
//			Line2D& line1 = vec_left_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_left_tar_line2d[line_serial_number2];
//			Line2D& line3 = vec_left_tar_line2d[line_serial_number3];
//			Point2f cross_point1 = conputeIntersectionPoint(line1, line2);
//			Point2f cross_point2 = conputeIntersectionPoint(line1, line3);
//			Point2f cross_point3 = conputeIntersectionPoint(line2, line3);
//			Point2f cross_point = (cross_point1 + cross_point2 + cross_point3) / 3.0;
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			if (point3.is_start_point == true)
//				line3.start_point = cross_point;
//			else
//				line3.end_point = cross_point;
//
//			cv::circle(left_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//	}
//	for (const auto& near_points : vec_right_near_points_group) // ��
//	{
//		if (near_points.size() == 1)
//			continue;
//		if (near_points.size() == 2)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			LinePoint2d point1 = vec_right_line_points[point_index1];
//			LinePoint2d point2 = vec_right_line_points[point_index2];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			Line2D& line1 = vec_right_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_right_tar_line2d[line_serial_number2];
//			Point2f cross_point = conputeIntersectionPoint(line1, line2);
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			cv::circle(right_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//		if (near_points.size() == 3)
//		{
//			int point_index1 = near_points[0];
//			int point_index2 = near_points[1];
//			int point_index3 = near_points[2];
//			LinePoint2d point1 = vec_right_line_points[point_index1];
//			LinePoint2d point2 = vec_right_line_points[point_index2];
//			LinePoint2d point3 = vec_right_line_points[point_index3];
//			int line_serial_number1 = point1.line_serial_number;
//			int line_serial_number2 = point2.line_serial_number;
//			int line_serial_number3 = point3.line_serial_number;
//			Line2D& line1 = vec_right_tar_line2d[line_serial_number1];
//			Line2D& line2 = vec_right_tar_line2d[line_serial_number2];
//			Line2D& line3 = vec_right_tar_line2d[line_serial_number3];
//			Point2f cross_point1 = conputeIntersectionPoint(line1, line2);
//			Point2f cross_point2 = conputeIntersectionPoint(line1, line3);
//			Point2f cross_point3 = conputeIntersectionPoint(line2, line3);
//			Point2f cross_point = (cross_point1 + cross_point2 + cross_point3) / 3.0;
//
//			if (point1.is_start_point == true)
//				line1.start_point = cross_point;
//			else
//				line1.end_point = cross_point;
//			if (point2.is_start_point == true)
//				line2.start_point = cross_point;
//			else
//				line2.end_point = cross_point;
//
//			if (point3.is_start_point == true)
//				line3.start_point = cross_point;
//			else
//				line3.end_point = cross_point;
//
//			cv::circle(right_cross_point_image, cross_point, 5, Scalar(0, 0, 255));
//		}
//	}
//
//	//// ��ʾֱ�߽���
//	//imshow("left_cross_point_image", left_cross_point_image);
//	//imshow("right_cross_point_image", right_cross_point_image);
//	//waitKey(0);
//
//	//// ��ʾ����������Ϊ�˵��ֱ��
//	//Mat left_line_img = left_tar_img.clone();
//	//Mat right_line_img = right_tar_img.clone();
//	//cvtColor(left_line_img, left_line_img, COLOR_GRAY2BGR);
//	//cvtColor(right_line_img, right_line_img, COLOR_GRAY2BGR);
//	//for (const auto& line2d : vec_left_tar_line2d)
//	//{
//	//	cv::line(left_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(left_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_tar_line2d)
//	//{
//	//	cv::line(right_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 1);
//	//	cv::putText(right_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//imshow("left_lines", left_line_img);
//	//imshow("right_lines", right_line_img);
//	//waitKey(0);
//
//
//	// ֱ�߾ֲ�����ӳ�䵽ȫ�����꣬���¼���б��
//	vector<Line2D> vec_left_line2d;
//	vector<Line2D> vec_right_line2d;
//	for (const auto& line2d : vec_left_tar_line2d) // ��
//	{
//		static int index = 0;
//		float x1 = line2d.start_point.x + rect_left_tar.x;
//		float y1 = line2d.start_point.y + rect_left_tar.y;
//		float x2 = line2d.end_point.x + rect_left_tar.x;
//		float y2 = line2d.end_point.y + rect_left_tar.y;
//		Vec4f line4f(x1, y1, x2, y2);
//		vec_left_line2d.push_back(Line2D(index, line4f));
//		++index;
//	}
//	for (const auto& line2d : vec_right_tar_line2d) // ��
//	{
//		static int index = 0;
//		float x1 = line2d.start_point.x + rect_right_tar.x;
//		float y1 = line2d.start_point.y + rect_right_tar.y;
//		float x2 = line2d.end_point.x + rect_right_tar.x;
//		float y2 = line2d.end_point.y + rect_right_tar.y;
//		Vec4f line4f(x1, y1, x2, y2);
//		vec_right_line2d.push_back(Line2D(index, line4f));
//		++index;
//	}
//
//	//// ���ƴ���õ�ȫ��ֱ��
//	//Mat left_global_line_img = rleftimg.clone();
//	//Mat right_global_line_img = rrightimg.clone();
//
//	//for (const auto& line2d : vec_left_line2d)
//	//{
//	//	cv::line(left_global_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 2);
//	//	cv::putText(left_global_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//	//for (const auto& line2d : vec_right_line2d)
//	//{
//	//	cv::line(right_global_line_img, line2d.start_point, line2d.end_point, Scalar(0, 0, 255), 2);
//	//	cv::putText(right_global_line_img, std::to_string(line2d.serial_number), line2d.mid_point, FONT_HERSHEY_SIMPLEX, 0.8, Scalar(255, 0, 0));
//	//}
//
//	//Mat left_right_global_line_img; // ƴ������ͼ��
//	//left_right_global_line_img.create(img_Height, img_Width * 2, CV_8UC3);
//	//hconcat(left_global_line_img, right_global_line_img, left_right_global_line_img);
//	//imshow("left_global_line_img", left_global_line_img);
//	//imshow("right_global_line_img", right_global_line_img);
//	//imshow("left_right_global_line_img", left_right_global_line_img);
//	//waitKey(0);
//
//
//	// ֱ��ƥ��
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_reconstruct(new pcl::PointCloud<pcl::PointXYZ>()); // ֱ�ߵ���
//	std::vector<std::vector<int>> vec_match_index_group; // ƥ������
//	lineMatch(vec_left_line2d, vec_right_line2d, vec_match_index_group);
//	
//	// ֱ���ؽ�
//	vector<Line3D> vec_line3d;
//	for (const auto vec_match_index : vec_match_index_group)
//	{
//		if (vec_match_index.size() != 2)
//			continue;
//
//		int left_index = vec_match_index[0];
//		int right_index = vec_match_index[1];
//		const Line2D& left_line = vec_left_line2d[left_index];
//		const Line2D& right_line = vec_right_line2d[right_index];
//
//		// p1��㣬p2�յ�
//		const Point2f& left_p1 = left_line.start_point;
//		const Point2f& left_p2 = left_line.end_point;
//		const Point2f& right_p1 = right_line.start_point;
//		const Point2f& right_p2 = right_line.end_point;
//
//		// ����ռ�ֱ�������˵�
//		const float f = 622.812; // ���࣬���㣬���߾���
//		const float cx = 654.193;
//		const float cy = 344.569;
//		const float Tx = 60.14;
//
//		float p1_disp = left_p1.x - right_p1.x;
//
//		float p1_z_w = f * Tx / p1_disp;
//		float p1_x_w = (left_p1.x - cx)*p1_z_w / f;
//		float p1_y_w = (left_p1.y - cy)*p1_z_w / f;
//
//		float p2_disp = left_p2.x - right_p2.x;
//		float p2_z_w = f * Tx / p2_disp;
//		float p2_x_w = (left_p2.x - cx)*p2_z_w / f;
//		float p2_y_w = (left_p2.y - cy)*p2_z_w / f;
//
//		// push�ؽ��Ŀռ�ֱ��
//		Line3D line_3d(Point3f(p1_x_w, p1_y_w, p1_z_w), Point3f(p2_x_w, p2_y_w, p2_z_w));
//		vec_line3d.push_back(line_3d);
//
//		// push �ռ�ֱ�߶˵�
//		pcl::PointXYZ p1_w(p1_x_w, p1_y_w, p1_z_w);
//		pcl::PointXYZ p2_w(p2_x_w, p2_y_w, p2_z_w);
//		cloud_reconstruct->points.push_back(p1_w);
//		cloud_reconstruct->points.push_back(p2_w);
//	}
//
//	// ȫ��ֱ�߲����ɵ���
//	float sample_dis = 1.0; // ֱ�ߴ�Ͼ���
//
//	for (const auto& line_3d : vec_line3d)
//	{
//		// ֱ�ߵȾ�����
//		float x0 = line_3d.x0_vec.x;
//		float y0 = line_3d.x0_vec.y;
//		float z0 = line_3d.x0_vec.z;
//
//		float v1 = line_3d.normal_vec.x;
//		float v2 = line_3d.normal_vec.y;
//		float v3 = line_3d.normal_vec.z;
//
//		float v_square = line_3d.normal_vec.dot(line_3d.normal_vec);
//		float deta_t = sample_dis / std::sqrt(v_square);
//
//		for (float t = line_3d.t_start; t < line_3d.t_end; t += deta_t)
//		{
//			pcl::PointXYZ point_3d;
//			point_3d.x = x0 + v1 * t;
//			point_3d.y = y0 + v2 * t;
//			point_3d.z = z0 + v3 * t;
//			cloud_reconstruct->points.push_back(point_3d);
//		}
//	}
//
//	// ������׼
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src(new pcl::PointCloud<pcl::PointXYZ>());
//	createCuboidWireFramePointCloud(cloud_src);
//	// ������Ʊ��淨��
//	clock_t fpfh_start_time = clock();
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_src; // Դ
//	ne_src.setInputCloud(cloud_src);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_src.setSearchMethod(tree_src);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_src.setRadiusSearch(2.0);
//	ne_src.setNumberOfThreads(8);
//	ne_src.compute(*cloud_src_normals);
//
//	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> ne_tar; // �ؽ�
//	ne_tar.setInputCloud(cloud_reconstruct);
//	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_reconstruct(new pcl::search::KdTree< pcl::PointXYZ>());
//	ne_tar.setSearchMethod(tree_reconstruct);
//	pcl::PointCloud<pcl::Normal>::Ptr cloud_reconstruct_normals(new pcl::PointCloud< pcl::Normal>);
//	ne_tar.setRadiusSearch(2.0);
//	ne_src.setNumberOfThreads(8);
//	ne_tar.compute(*cloud_reconstruct_normals);
//
//	//����FPFH
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_src; // Դ
//	fpfh_src.setInputCloud(cloud_src);
//	fpfh_src.setInputNormals(cloud_src_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_src_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_src.setSearchMethod(tree_src_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_src(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_src.setRadiusSearch(2.0);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_src.compute(*fpfhs_src);
//	std::cout << "compute *cloud_src fpfh" << endl;
//
//	pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_tar; // �ؽ�
//	fpfh_tar.setInputCloud(cloud_reconstruct);
//	fpfh_tar.setInputNormals(cloud_reconstruct_normals);
//	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_reconstruct_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
//	fpfh_tar.setSearchMethod(tree_reconstruct_fpfh);
//	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_reconstruct(new pcl::PointCloud<pcl::FPFHSignature33>());
//	fpfh_tar.setRadiusSearch(2.0);
//	fpfh_src.setNumberOfThreads(8);
//	fpfh_tar.compute(*fpfhs_reconstruct);
//
//	std::cout << "compute *cloud_tar fpfh" << endl;
//	clock_t fpfh_end_time = clock();
//	cout << "fpfh time: " << (double)(fpfh_end_time - fpfh_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//SAC��׼
//	clock_t sac_start_time = clock();
//
//	pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> scia;
//	scia.setInputSource(cloud_src);
//	scia.setInputTarget(cloud_reconstruct);
//	scia.setSourceFeatures(fpfhs_src);
//	scia.setTargetFeatures(fpfhs_reconstruct);
//
//	//scia.setMinSampleDistance(1);
//	//scia.setNumberOfSamples(2);
//	//scia.setCorrespondenceRandomness(20);
//
//	pcl::PointCloud<pcl::PointXYZ>::Ptr sac_result(new pcl::PointCloud<pcl::PointXYZ>);
//	scia.align(*sac_result);
//	Eigen::Matrix4f sac_trans = scia.getFinalTransformation();
//	clock_t sac_end_time = clock();
//	cout << "sac time: " << (double)(sac_end_time - sac_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//
//	//icp��׼
//	clock_t icp_start_time = clock();
//	pcl::PointCloud<pcl::PointXYZ>::Ptr icp_result(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
//	icp.setInputSource(cloud_src);
//	icp.setInputTarget(cloud_reconstruct);
//	// ��������Ӧ����� (e.g., correspondences with higher distances will be ignored)
//	icp.setMaxCorrespondenceDistance(40);
//	// ����������
//	icp.setMaximumIterations(100);
//	// ���α仯����֮��Ĳ�ֵ
//	icp.setTransformationEpsilon(1e-10);
//	// �������
//	icp.setEuclideanFitnessEpsilon(0.2);
//	icp.align(*icp_result, sac_trans);
//	clock_t icp_end_time = clock();
//	cout << "icp time: " << (double)(icp_end_time - icp_start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
//	// ��ȡ���յı任����
//	Eigen::Matrix4f final_trans_pcl = icp.getFinalTransformation();
//	cv::Mat_<float> final_trans_cv(4, 4);
//	for (size_t i = 0; i < 4; i++)
//	{
//		for (size_t j = 0; j < 4; j++)
//		{
//			final_trans_cv(i, j) = final_trans_pcl(i, j);
//		}
//	}
//
//	// ��������ϵ������
//	vector<cv::Point3f> vec_axis_3d;
//	vector<cv::Point2f> vec_axis_2d;
//	vec_axis_3d.push_back(Point3f(10.0, 0.0, 0.0));
//	vec_axis_3d.push_back(Point3f(0.0, 10.0, 0.0));
//	vec_axis_3d.push_back(Point3f(0.0, 0.0, 10.0));
//	vec_axis_3d.push_back(Point3f(0.0, 0.0, 0.0));
//
//	// ��������ϵת��������
//	Eigen::Vector3f final_angle;
//	Eigen::Vector3f final_translate;
//	matrix2Angle(final_trans_pcl, final_angle);
//	matrix2Translate(final_trans_pcl, final_translate);
//	cv::Mat rVec(3, 1, cv::DataType<float>::type); // Rotation vector
//	rVec.at<float>(0) = final_angle.x();
//	rVec.at<float>(1) = final_angle.y();
//	rVec.at<float>(2) = final_angle.z();
//
//	cv::Mat tVec(3, 1, cv::DataType<float>::type); // Translation vector
//	tVec.at<float>(0) = final_translate.x();
//	tVec.at<float>(1) = final_translate.y();
//	tVec.at<float>(2) = final_translate.z();
//
//	cout << "λ�ƣ�" << tVec << endl;
//	cout << "��̬��" << rVec << endl;
//
//	cv::projectPoints(vec_axis_3d, rVec, tVec, Pl.colRange(0, 3), Mat(), vec_axis_2d);
//
//	// ������ͶӰ���������
//	Mat reproject_axis_img = rleftimg.clone();
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[0], Scalar(0, 0, 255), 2);
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[1], Scalar(0, 255, 0), 2);
//	cv::line(reproject_axis_img, vec_axis_2d[3], vec_axis_2d[2], Scalar(255, 0, 0), 2);
//	imshow("reproject_axis_img", reproject_axis_img);
//	waitKey(0);
//
//	// ��ʾ����
//	// ��ɫ��Դ����
//	// ��ɫ��Ŀ�����
//	// ��ɫ����׼����
//	pcl::visualization::PCLVisualizer viewer("registration Viewer");
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h(cloud_src, 0, 255, 0); // ģ��RGB
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> tgt_h(cloud_reconstruct, 255, 0, 0); // �ؽ���
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> final_h(icp_result, 0, 0, 255); // ��׼���ģ��
//	viewer.setBackgroundColor(1.0, 1.0, 1.0);
//
//	//viewer.addPointCloud(cloud_src, src_h, "tar_points");
//	viewer.addPointCloud(cloud_reconstruct, tgt_h, "src_points");
//	viewer.addPointCloud(icp_result, final_h, "icp_points");
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "src_points"); // ���õ��ƴ�С���ؽ���
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "tar_points"); // ģ��
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "icp_points"); // ��׼���ģ��
//	viewer.addCoordinateSystem(20.0); // �������ָʾ�ᣬ���ű���-------��ɫX����ɫY����ɫZ
//
//	//// �����ؽ�����
//	//cloud_reconstruct->width = 1;
//	//cloud_reconstruct->height = cloud_reconstruct->points.size();
//	//pcl::io::savePCDFileASCII("./cloud_reconstruct_src/c_reconstruct.pcd", *cloud_reconstruct);
//	//cout << "cloud_reconstruct is saved" << endl;
//
//	//// ����Դ����
//	//cloud_src->width = 1;
//	//cloud_src->height = cloud_src->points.size();
//	//pcl::io::savePCDFileASCII("./cloud_reconstruct_src/c_src.pcd", *cloud_src);
//	//cout << "cloud_src is saved" << endl;
//
//	//// ������׼����
//	//icp_result->width = 1;
//	//icp_result->height = icp_result->points.size();
//	//pcl::io::savePCDFileASCII("./cloud_reconstruct_src/c_icp_result.pcd", *icp_result);
//	//cout << "icp_result is saved" << endl;
//
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(100);
//		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	}
//
//	return 0;
//}
//
////����תƽ�ƾ��������ת�Ƕ�
//void matrix2Angle(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_angle)
//{
//	double ax, ay, az;
//	if (result_trans(2, 0) == 1 || result_trans(2, 0) == -1)
//	{
//		az = 0;
//		double dlta;
//		dlta = atan2(result_trans(0, 1), result_trans(0, 2));
//		if (result_trans(2, 0) == -1)
//		{
//			ay = M_PI / 2;
//			ax = az + dlta;
//		}
//		else
//		{
//			ay = -M_PI / 2;
//			ax = -az + dlta;
//		}
//	}
//	else
//	{
//		ay = -asin(result_trans(2, 0));
//		ax = atan2(result_trans(2, 1) / cos(ay), result_trans(2, 2) / cos(ay));
//		az = atan2(result_trans(1, 0) / cos(ay), result_trans(0, 0) / cos(ay));
//	}
//	result_angle << ax, ay, az;
//}
//
////����תƽ�ƾ������ƽ�ƾ���
//void matrix2Translate(Eigen::Matrix4f &result_trans, Eigen::Vector3f &result_translate)
//{
//	result_translate << result_trans(0, 3), result_trans(1, 3), result_trans(2, 3);
//}
//
//// �����������߿����
//void createCuboidWireFramePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_wire_frame)
//{
//
//	int step = 1;
//
//	// ���ɳ�
//	for (int i = 0; i <= 50; i += step)
//	{
//		pcl::PointXYZ p1;
//		p1.x = i; p1.y = 0; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = i; p2.y = 40; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = i; p3.y = 0; p3.z = 30;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = i; p4.y = 40; p4.z = 30;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//
//	// ���ɿ�
//	for (int i = 0; i <= 40; i += step)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 0; p1.y = i; p1.z = 0;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 50; p2.y = i; p2.z = 0;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 0; p3.y = i; p3.z = 30;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 50; p4.y = i; p4.z = 30;
//		cloud_wire_frame->points.push_back(p4);
//	}
//
//	// ���ɸ�
//	for (int i = 0; i <= 30; i += step)
//	{
//		pcl::PointXYZ p1;
//		p1.x = 0; p1.y = 0; p1.z = i;
//		cloud_wire_frame->points.push_back(p1);
//
//		pcl::PointXYZ p2;
//		p2.x = 50; p2.y = 0; p2.z = i;
//		cloud_wire_frame->points.push_back(p2);
//
//		pcl::PointXYZ p3;
//		p3.x = 0; p3.y = 40; p3.z = i;
//		cloud_wire_frame->points.push_back(p3);
//
//		pcl::PointXYZ p4;
//		p4.x = 50; p4.y = 40; p4.z = i;
//		cloud_wire_frame->points.push_back(p4);
//	}
//}


///* ��ȡ����ʵ���pcd�ļ� �������յ�ʵ���� */
//#include <iostream>
//
//#include <opencv2/opencv.hpp>
//#include <opencv2/ximgproc.hpp>
//#include <opencv2/xfeatures2d.hpp>
//
//#include <pcl/point_cloud.h>
//#include <pcl/io/pcd_io.h>
//#include <pcl/point_types.h>
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/common/transforms.h>
//
//int user_data;
//
//// ͼƬ��С
//const int img_Width = 1280;
//const int img_Height = 720;
//
////У����ת����R��ͶӰ����P ��ͶӰ����Q
//cv::Mat Rl, Rr, Pl, Pr, Q;
//
//int main(int argc, char** argv)
//{
//
//	cv::Mat left_and_right = cv::imread("./left_right_img/IMG_0022.jpg");
//	cv::Rect rect_left(0, 0, img_Width, img_Height);
//	cv::Rect rect_right(img_Width, 0, img_Width, img_Height);
//
//	cv::Mat left = left_and_right(rect_left);
//	cv::Mat right = left_and_right(rect_right);
//
//	//imshow("src_left", left);
//	//imshow("src_right", right);
//	//waitKey(0);
//
//	cv::Mat cameraMatrix[2], distCoeffs[2];
//	cv::Mat Q;
//
//	// ����ڲ�
//	cameraMatrix[0] = cv::Mat((cv::Mat_<double>(3, 3) << 731.5022, 0.0, 0.0,
//		0.0, 732.0644, 0.0,
//		674.6773, 357.5894, 1.0)).t();
//
//	cameraMatrix[1] = cv::Mat((cv::Mat_<double>(3, 3) << 727.3993, 0.0, 0.0,
//		0.0, 727.6805, 0.0,
//		634.0939, 339.6110, 1.0)).t();
//
//	// ����ϵ��
//	distCoeffs[0] = (cv::Mat_<double>(4, 1) << 0.1070, -0.1208, 0.0, -0.0004);
//	distCoeffs[1] = (cv::Mat_<double>(4, 1) << 0.1145, -0.1517, 0.0007, -0.0004);
//
//	// ������λ��
//	cv::Mat R = cv::Mat((cv::Mat_<double>(3, 3) << 1.0, 4.7372e-04, 4.8307e-04,
//		-4.7146e-04, 1.0, -0.0047,
//		-0.0005, 0.0047, 1.0)).t();
//	cv::Mat T = (cv::Mat_<double>(3, 1) << -60.1450, -0.0701, -0.0956);
//
//	// ͼ��У��֮�����Ч����
//	cv::Rect validRoi[2];
//
//	// bougust����У��
//	cv::stereoRectify(cameraMatrix[0], distCoeffs[0],
//		cameraMatrix[1], distCoeffs[1],
//		cv::Size(img_Width, img_Height), R, T, Rl, Rr, Pl, Pr, Q,
//		cv::CALIB_ZERO_DISPARITY, 1, cv::Size(img_Width, img_Height), &validRoi[0], &validRoi[1]); //validPixROI1  validPixROI2 - ��ѡ�����������Rect�����ݡ����ڲ����������ض���Ч
//
//	//��ʾ���������ͶӰ����
//	cout << "camL" << cameraMatrix[0] << endl;
//	cout << "camR" << cameraMatrix[1] << endl;
//	cout << "Pro1" << Pl << endl;
//	cout << "Pro2" << Pr << endl;
//	cout << "Q" << Q << endl;
//
//	// ����������ͼ��У������ӳ���
//	cv::Mat rmap[2][2];
//	cv::initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], Rl, Pl, cv::Size(img_Width, img_Height), CV_16SC2, rmap[0][0], rmap[0][1]);
//	cv::initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], Rr, Pr, cv::Size(img_Width, img_Height), CV_16SC2, rmap[1][0], rmap[1][1]);
//
//	// ��Դͼ�������ӳ�䣬ӳ�����ΪУ�����ͼ��
//	cv::Mat rleftimg, rrightimg, rcolor;
//	cv::remap(left, rleftimg, rmap[0][0], rmap[0][1], cv::INTER_LINEAR);
//	cv::remap(right, rrightimg, rmap[1][0], rmap[1][1], cv::INTER_LINEAR);
//
//	////// ����У����ͼƬ
//	////imwrite("rleft.jpg", rleftimg);
//	////imwrite("rright.jpg", rrightimg);
//
//	//// �������߽��жԱ�
//	//cv::Mat canvas;
//	//canvas.create(img_Height, img_Width * 2, CV_8UC3);
//	//hconcat(rleftimg, rrightimg, canvas);
//	//for (int j = 0; j < canvas.rows; j += 16)
//	//	cv::line(canvas, cv::Point(0, j), cv::Point(canvas.cols, j), cv::Scalar(0, 255, 0), 1, 8);
//	//cv::line(canvas, cv::Point(1280, 0), cv::Point(1280, 720), cv::Scalar(255, 255, 0), 1, 8);
//
//	//// ��ʾУ�����
//	//cv::resize(canvas, canvas, cv::Size(), 0.5, 0.5);
//	//cv::namedWindow("rectifyResult", cv::WINDOW_AUTOSIZE);
//	//cv::imshow("rectifyResult", canvas);
//	//cv::waitKey(0);
//
//	// // ��ʾУ�����ƴ������ͼ��
//	//cv::Mat left_right_global_line_img;
//	//left_right_global_line_img.create(img_Height, img_Width * 2, CV_8UC3);
//	//cv::hconcat(rleftimg, rrightimg, left_right_global_line_img);
//	//cv::imshow("left_right_global_line_img", left_right_global_line_img);
//	//cv::waitKey(0);
//
//	// ��ȡ�����ļ�
//	pcl::PointCloud<pcl::PointXYZ>::Ptr c_src_cloud(new pcl::PointCloud<pcl::PointXYZ>); // ������
//	pcl::PointCloud<pcl::PointXYZ>::Ptr c_reconstruct_cloud(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::io::loadPCDFile<pcl::PointXYZ>("./cloud_reconstruct_src/c_src.pcd", *c_src_cloud);
//	pcl::io::loadPCDFile<pcl::PointXYZ>("./cloud_reconstruct_src/c_reconstruct.pcd", *c_reconstruct_cloud);
//
//	pcl::PointCloud<pcl::PointXYZ>::Ptr t_src_cloud(new pcl::PointCloud<pcl::PointXYZ>); // T
//	pcl::PointCloud<pcl::PointXYZ>::Ptr t_reconstruct_cloud(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::io::loadPCDFile<pcl::PointXYZ>("./cloud_reconstruct_src/t_src.pcd", *t_src_cloud);
//	pcl::io::loadPCDFile<pcl::PointXYZ>("./cloud_reconstruct_src/t_reconstruct.pcd", *t_reconstruct_cloud);
//
//	pcl::PointCloud<pcl::PointXYZ>::Ptr l_src_cloud(new pcl::PointCloud<pcl::PointXYZ>); // L
//	pcl::PointCloud<pcl::PointXYZ>::Ptr l_reconstruct_cloud(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::io::loadPCDFile<pcl::PointXYZ>("./cloud_reconstruct_src/l_src.pcd", *l_src_cloud);
//	pcl::io::loadPCDFile<pcl::PointXYZ>("./cloud_reconstruct_src/l_reconstruct.pcd", *l_reconstruct_cloud);
//
//
//	// ������ͶӰ���������
//	cv::Mat reproject_axis_img = rleftimg.clone();
//	cv::imwrite("reproject_axis_img.jpg", reproject_axis_img);
//	std::vector<cv::Point3f> vec_axis_3d;
//	vec_axis_3d.push_back(cv::Point3f(20.0, 0.0, 0.0));
//	vec_axis_3d.push_back(cv::Point3f(0.0, 20.0, 0.0));
//	vec_axis_3d.push_back(cv::Point3f(0.0, 0.0, 20.0));
//	vec_axis_3d.push_back(cv::Point3f(0.0, 0.0, 0.0));
//
//	std::vector<cv::Point2f> c_vec_axis_2d; // ������
//	cv::Mat c_rVec(3, 1, cv::DataType<float>::type);
//	c_rVec.at<float>(0) = 2.5913;
//	c_rVec.at<float>(1) = 0.2678;
//	c_rVec.at<float>(2) = 0.3763;
//	cout << "������" << c_rVec * 180.0 / CV_PI << endl;
//	cv::Mat c_tVec(3, 1, cv::DataType<float>::type);
//	c_tVec.at<float>(0) = -85.8915;
//	c_tVec.at<float>(1) = -1.1579;
//	c_tVec.at<float>(2) = 199.0901;
//	cv::projectPoints(vec_axis_3d, c_rVec, c_tVec, Pl.colRange(0, 3), cv::Mat(), c_vec_axis_2d);
//	cv::line(reproject_axis_img, cv::Point2f(502, 399), cv::Point2f(444, 358), cv::Scalar(0, 0, 255), 2);
//	cv::line(reproject_axis_img, cv::Point2f(502, 399), cv::Point2f(553, 357), cv::Scalar(0, 255, 0), 2);
//	cv::line(reproject_axis_img, cv::Point2f(502, 399), cv::Point2f(487, 349), cv::Scalar(255, 0, 0), 2);
//
//	std::vector<cv::Point2f> t_vec_axis_2d; // T
//	cv::Mat t_rVec(3, 1, cv::DataType<float>::type);
//	t_rVec.at<float>(0) = 2.4176;
//	t_rVec.at<float>(1) = 0.0204;
//	t_rVec.at<float>(2) = -0.025;
//	cout << "Tƽ��Ƭ" << t_rVec * 180.0 / CV_PI << endl;
//	cv::Mat t_tVec(3, 1, cv::DataType<float>::type);
//	t_tVec.at<float>(0) = 19.0347;
//	t_tVec.at<float>(1) = 55.5203;
//	t_tVec.at<float>(2) = 144.8861;
//	cv::projectPoints(vec_axis_3d, t_rVec, t_tVec, Pl.colRange(0, 3), cv::Mat(), t_vec_axis_2d);
//	cv::line(reproject_axis_img, t_vec_axis_2d[3], t_vec_axis_2d[0], cv::Scalar(0, 0, 255), 2);
//	cv::line(reproject_axis_img, t_vec_axis_2d[3], t_vec_axis_2d[1], cv::Scalar(0, 255, 0), 2);
//	cv::line(reproject_axis_img, t_vec_axis_2d[3], t_vec_axis_2d[2], cv::Scalar(255, 0, 0), 2);
//
//	std::vector<cv::Point2f> l_vec_axis_2d; // L
//	cv::Mat l_rVec(3, 1, cv::DataType<float>::type);
//	l_rVec.at<float>(0) = 2.4316;
//	l_rVec.at<float>(1) = -0.1272;
//	l_rVec.at<float>(2) = -0.2181;
//	cout << "L�ͽ���" << l_rVec * 180.0 / CV_PI << endl;
//	cv::Mat l_tVec(3, 1, cv::DataType<float>::type);
//	l_tVec.at<float>(0) = 85.3124;
//	l_tVec.at<float>(1) = 19.6473;
//	l_tVec.at<float>(2) = 177.3346;
//	cv::projectPoints(vec_axis_3d, l_rVec, l_tVec, Pl.colRange(0, 3), cv::Mat(), l_vec_axis_2d);
//	cv::line(reproject_axis_img, l_vec_axis_2d[3], l_vec_axis_2d[0], cv::Scalar(0, 0, 255), 2);
//	cv::line(reproject_axis_img, l_vec_axis_2d[3], l_vec_axis_2d[1], cv::Scalar(0, 255, 0), 2);
//	cv::line(reproject_axis_img, l_vec_axis_2d[3], l_vec_axis_2d[2], cv::Scalar(255, 0, 0), 2);
//
//	imshow("reproject_axis_img", reproject_axis_img);
//	cv::waitKey(0);
//
//
//
//	// �任����
//	Eigen::Affine3f c_transform = Eigen::Affine3f::Identity();  // ������
//	c_transform.translation() << -85.8915, -1.1579, 199.0901; float theta_z = 0.3763; float theta_x = 2.5913; float theta_y = 0.2678;
//	c_transform.rotate(Eigen::AngleAxisf(theta_z, Eigen::Vector3f::UnitZ()));
//	c_transform.rotate(Eigen::AngleAxisf(theta_y, Eigen::Vector3f::UnitY()));
//	c_transform.rotate(Eigen::AngleAxisf(theta_x, Eigen::Vector3f::UnitX()));
//	pcl::PointCloud<pcl::PointXYZ>::Ptr c_transformed_cloud(new pcl::PointCloud<pcl::PointXYZ>());
//	pcl::transformPointCloud(*c_src_cloud, *c_src_cloud, c_transform);
//
//	Eigen::Affine3f t_transform = Eigen::Affine3f::Identity();  // T
//	t_transform.translation() << 19.0347, 55.5203, 144.8861;
//	t_transform.rotate(Eigen::AngleAxisf(-0.025, Eigen::Vector3f::UnitZ()));
//	t_transform.rotate(Eigen::AngleAxisf(0.0204, Eigen::Vector3f::UnitY()));
//	t_transform.rotate(Eigen::AngleAxisf(2.4176, Eigen::Vector3f::UnitX()));
//	pcl::PointCloud<pcl::PointXYZ>::Ptr t_transformed_cloud(new pcl::PointCloud<pcl::PointXYZ>());
//	pcl::transformPointCloud(*t_src_cloud, *t_src_cloud, t_transform);
//
//	Eigen::Affine3f l_transform = Eigen::Affine3f::Identity();  // L
//	l_transform.translation() << 85.3124, 19.6473, 177.3346;
//	l_transform.rotate(Eigen::AngleAxisf(-0.2181, Eigen::Vector3f::UnitZ()));
//	l_transform.rotate(Eigen::AngleAxisf(-0.1272, Eigen::Vector3f::UnitY()));
//	l_transform.rotate(Eigen::AngleAxisf(2.4316, Eigen::Vector3f::UnitX()));
//	pcl::PointCloud<pcl::PointXYZ>::Ptr l_transformed_cloud(new pcl::PointCloud<pcl::PointXYZ>());
//	pcl::transformPointCloud(*l_src_cloud, *l_src_cloud, l_transform);
//
//
//	// ��ʾ����
//	pcl::visualization::PCLVisualizer viewer("registration Viewer");
//	viewer.setBackgroundColor(1.0, 1.0, 1.0);
//	viewer.addCoordinateSystem(20.0); // �������ָʾ�ᣬ���ű���-------��ɫX����ɫY����ɫZ
//
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> c_reconstruct_cloud_h(c_reconstruct_cloud, 255, 0, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> c_src_cloud_h(c_src_cloud, 0, 0, 255);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> t_reconstruct_cloud_h(t_reconstruct_cloud, 255, 0, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> t_src_cloud_h(t_src_cloud, 0, 0, 255);	
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> l_reconstruct_cloud_h(l_reconstruct_cloud, 255, 0, 0);
//	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> l_src_cloud_h(l_src_cloud, 0, 0, 255);
//
//	viewer.addPointCloud(c_reconstruct_cloud, c_reconstruct_cloud_h, "c_reconstruct_points");
//	viewer.addPointCloud(c_src_cloud, c_src_cloud_h, "c_src_points");
//	viewer.addPointCloud(t_reconstruct_cloud, t_reconstruct_cloud_h, "t_reconstruct_points");
//	viewer.addPointCloud(t_src_cloud, t_src_cloud_h, "t_src_points");
//	viewer.addPointCloud(l_reconstruct_cloud, l_reconstruct_cloud_h, "l_reconstruct_points");
//	viewer.addPointCloud(l_src_cloud, l_src_cloud_h, "l_src_points");
//
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "c_reconstruct_points"); // ���õ��ƴ�С
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "c_src_points");
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "t_reconstruct_points");
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "t_src_points"); 
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "l_reconstruct_points");
//	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "l_src_points");
//
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(100);
//		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
//	}
//
//	return 0;
//}