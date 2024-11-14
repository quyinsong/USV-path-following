# USV-path-following
USV路径跟踪LOS控制算法仿真

实验一(main1)： 直线路径跟踪，LOS导航算法 PID控制
1、USV数学模型

参考文献；Error_Constrained LOS Path Following of a Surface Vessel With Actuator Saturation and Faults\
USV 数学模型矩阵表达式如下：
	![image](https://github.com/user-attachments/assets/a7d36b02-a73f-4a7d-aeb5-cc653a50f9f7) （1）
其中：

![image](https://github.com/user-attachments/assets/869b428b-ef0c-4786-8c8e-01cc3d6d5957)

USV仿真参数：船长：1.255m 船宽： 0.29m
			
![image](https://github.com/user-attachments/assets/e0427c56-9cca-4c75-92df-1bdd3d0b4823)


可写成如下表达式：

![image](https://github.com/user-attachments/assets/d54f75e2-df43-4e71-b409-771f4e006779)

![image](https://github.com/user-attachments/assets/d3b10606-6ea1-4aee-8720-8f37d11f4cef)

其中：

![image](https://github.com/user-attachments/assets/361cf336-c65f-4cd6-ab53-2c7caf5870df)

2 LOS制导率

![image](https://github.com/user-attachments/assets/1a348499-3e8d-4ffb-8cd7-ebc9422d299b)

控制目标：

直线路径是由多个点连接形成，控制目标是使USV跟踪期望直线路径，使其横向偏差ye->0

![image](https://github.com/user-attachments/assets/ebb65024-0b69-4bec-9763-a84d79d75d00)

期望航向角

![image](https://github.com/user-attachments/assets/95689774-f7c1-458c-8609-a8d6969536e1)

期望艏向角

![image](https://github.com/user-attachments/assets/afa8ab90-b9e6-4692-8247-b2321178bf16)

证明：

![image](https://github.com/user-attachments/assets/21482c75-a2c2-4528-830f-9c0fc459ffda)

![image](https://github.com/user-attachments/assets/465d017c-7643-4757-b21a-ecc15c7e28d8)

利用李雅普诺夫定理可知系统稳定
证毕；
3 艏向控制率设计
由USV数学模型可得：

![image](https://github.com/user-attachments/assets/0f3417f9-5432-434a-a022-e57b07d55adb)

可以通过调节Kp和kd的值来改变系统的响应速度
4 仿真
LOS参数：
PD控制器参数：Kp = 4  Kd = 6
点集：point_database =[0 0; 40 40; 80 40; 90 20; 90 10; 80 0]';
仿真结果：

![image](https://github.com/user-attachments/assets/98e1b61e-5e54-41a2-b9ce-310f3dc31e39)
