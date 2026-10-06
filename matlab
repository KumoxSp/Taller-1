% 1. Lectura de datos
tabla = readtable("datos_motor.csv");
t = tabla{:,2};
u = tabla{:,3};
y = tabla{:,4};

% 2. Hallar K
delta_u = max(u);
[delta_y_max, lin_max] = max(y);

% EXTRAEMOS LOS VALORES MÁXIMOS Y MÍNIMOS EN VARIABLES
y_max_prom = mean([y(lin_max),y(lin_max-1),y(lin_max+1)]);
y_min_prom = mean(y(1:find(u,1)));

delta_y = y_max_prom - y_min_prom;
k = delta_y / delta_u;

% --- 3. Método Nichols (Tangente) ---
derivada = diff(y) ./ diff(t);
[max_pen, ubi] = max(derivada);
tangente = max_pen * (t - t(ubi)) + y(ubi);
theta_nic = t(ubi) - (y(ubi) / max_pen);
tao_nic = delta_y / max_pen;
s = tf('s');
G_nic = (k / (tao_nic*s + 1)) * exp(-theta_nic*s);
[y_sim_nic, t_sim_nic] = step(delta_u * G_nic, t); % Extraemos datos para graficar

% --- 4. Método Miller (63.2%, theta compartido) ---
por63 = delta_y * 0.632;
[~, er_ubi] = min(abs(y - por63));
tao_63 = t(er_ubi) - theta_nic;
G_63 = (k / (tao_63*s + 1)) * exp(-theta_nic*s);
[y_sim_63, t_sim_63] = step(delta_u * G_63, t);

% --- 5. Método Analítico (Broida: 28.3% y 63.2%) ---
por28 = 0.283 * delta_y; 
[~, idx_t1] = min(abs(y - por28)); 
[~, idx_t2] = min(abs(y - por63)); 
t1 = t(idx_t1); % Tiempo al 28.3%
t2 = t(idx_t2); % Tiempo al 63.2%
tao_ana = 1.5 * (t2 - t1);
theta_ana = t2 - tao_ana;
G_ana = (k / (tao_ana*s + 1)) * exp(-theta_ana*s);
[y_sim_ana, t_sim_ana] = step(delta_u * G_ana, t);

% --- 6. Graficar todos los resultados ---
figure;
hold on;

% Datos originales
plot(t, u, 'b', t, y, 'g', t, tangente, 'r', 'LineWidth', 1.5);

% Modelos simulados
plot(t_sim_nic, y_sim_nic, 'w-.', 'LineWidth', 2); % Nichols (Negro)
plot(t_sim_63, y_sim_63, 'm--', 'LineWidth', 2);   % Miller (Magenta)
plot(t_sim_ana, y_sim_ana, 'c:', 'LineWidth', 2);  % Analítico (Cian)

% NUEVO: Trazos horizontales de los límites máximo y mínimo
% Usamos un color gris [0.5 0.5 0.5] y línea punteada '--'
yline(y_max_prom, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5); 
% Al segundo le ponemos HandleVisibility='off' para que no duplique la leyenda
yline(y_min_prom, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 

% Formato de gráfica
xlabel('Tiempo (s)'); ylabel('Amplitud');
title('Comparación de Métodos de Identificación');

% Ajusté el ylim para que use delta_y_max en lugar de delta_y (evita cortes si la base no es 0)
ylim([0, delta_y_max * 1.2]); 

% ACTUALIZADA LA LEYENDA: Se añade 'Límites Max/Min' al final
legend('u(t)', 'y(t) Real', 'Tangente', 'Modelo Nichols', 'Modelo 63.2%', ...
    'Modelo Analítico', 'Límites Max/Min', 'Location', 'southeast');
grid on;
hold off;
